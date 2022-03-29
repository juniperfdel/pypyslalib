import argparse
import re
from pathlib import Path

import autopep8

file_args = []
inputted_files = set()
converted = set()


class Transformation:
	def __init__(self, f_test, py_replace):
		self.f_t = re.compile(f_test)
		self.py_r = py_replace
	
	def test(self, in_string):
		return self.f_t.search(in_string) is not None if in_string else False
	
	def replace(self, in_string):
		return self.f_t.sub(self.py_r, in_string)


class DependT(Transformation):
	def replace(self, in_string):
		rv = re.sub(self.f_t, self.py_r, in_string)
		depends_ = [
			x.replace("sla", "").replace("_", "").lower().strip()
			for x in rv.split(":")[1].strip().split(",")
		]
		depends_ = [
			x for x in depends_ if
			(len(x) > 0) and (x not in converted)
		]
		for x in depends_:
			if x not in converted:
				inputted_files.add(x)
		return rv


class CallT(Transformation):
	def replace(self, in_string):
		rv = re.sub(self.f_t, self.py_r, in_string)
		match = re.search(r"sla_(\w+)\(", rv)
		if not match:
			return rv
		rv = re.sub(
			r"sla_(\w+)\(", match.group(1).lower() + "(", rv
		)
		new_dep = match.group(1).lower()
		if new_dep not in converted:
			print("Adding Dependency ", new_dep)
			inputted_files.add(new_dep)
		return rv


class FuncT(Transformation):
	def replace(self, in_string):
		rv = re.sub(self.f_t, self.py_r, in_string)
		match = re.search(r"sla_(\w+)\((.+)\)", rv)
		fn_args = [x.strip() for x in match.group(2).split(",")]
		file_args.extend(fn_args)
		fn_args = [x.lower() for x in fn_args]
		fn_args_str = ", ".join(fn_args)
		rv = re.sub(
			r"sla_(\w+)\((.+)\)",
			match.group(1).lower() + "(" + fn_args_str + ")",
			rv
		)
		return rv


line_replace = {
	"function": FuncT(
		r"SUBROUTINE (\w+) ?\(([ A-Za-z0-9,]+)\)",
		r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):"
	),
	"f_function": FuncT(
		r".+FUNCTION (\w+) ?\(([ A-Za-z0-9,]+)\)",
		r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):"
	),
	"dependencies": DependT(
		r"^\s*\*+\s*Called:(.+)", r"# Depends:\g<1>"
	),
	"comment": Transformation(r"^\*+\s*(.*)", r"# \g<1>"),
	"call": CallT(r"CALL (\w+) ?\((.+)\)", r"cls.\g<1>(\g<2>)"),
	"do": Transformation(
		r"DO (\S+)=(\S+),(\S+)", r"for \g<1> in (\g<2>, \g<3>):"
	),
	"end": Transformation(r"END (DO|IF)", r""),
	"eq_sign": Transformation("=", "="),
	"elif": Transformation(r"ELSE IF", r"elif:"),
	"if": Transformation(r"IF(.+)\sTHEN", r"if\g<1>:"),
	"else": Transformation(r"ELSE", r"else:"),
	"sl_if": Transformation(r"IF \((.+)\)(.+)", r"if \g<1>: \g<2>"),
	"lt": Transformation(r" ?\.LT\. ?", r" < "),
	"gt": Transformation(r" ?\.GT\. ?", r" > "),
	"eq": Transformation(r" ?\.EQ\. ?", r" == "),
	"le": Transformation(r" ?\.LE\. ?", r" <= "),
	"ge": Transformation(r" ?\.GE\. ?", r" >= "),
	"neq": Transformation(r" ?\.NE\. ?", r" != "),
	"and": Transformation(r" ?\.AND\. ?", r" and "),
	"or": Transformation(r" ?\.OR\. ?", r" or "),
	"anint": Transformation(r"ANINT", r"np.rint"),
	"dble": Transformation(r"DBLE", ""),
	"sqrt": Transformation(r"SQRT", r"np.sqrt"),
	"sign": Transformation(r"SIGN", r"np.sign"),
	"nint": Transformation(r"NINT", r"np.rint"),
	"abs": Transformation(r"ABS", r"np.abs"),
	"cos": Transformation(r"COS\(", r"np.cos("),
	"sin": Transformation(r"SIN\(", r"np.sin("),
	"tan": Transformation(r"TAN\(", r"np.tan("),
	"atan2": Transformation(r"ATAN2\(", r"np.arctan2("),
	"mod": Transformation(r"MOD", r"np.mod"),
	"setter": Transformation(r"(\w+)\((\d+)\) = ", r"\g<1>[\g<2>] = "),
	"max": Transformation(r"MAX", r"np.maximum"),
	"min": Transformation(r"MIN", "np.minimum"),
	"end_file": Transformation("END", ""),
	"implicit_none": Transformation("IMPLICIT NONE", ""),
	"double_prec": Transformation(r"DOUBLE PRECISION (.+)", r""),
	"int": Transformation(r"INTEGER (.+)", r""),
	"double_constant": Transformation(
		r"(\d+)(?:D|E)(-?)\+?0*(\d+)", r"\g<1>e\g<2>\g<3>"
	),
	"file_return": Transformation(r"^sla_\w+ = (.+)", r"return \g<1>")
}

# (Current_line, Next_line)
indent_change = {
	"f_function": (0, 2),
	"function": (0, 2),
	"do": (0, 1),
	"if": (0, 1),
	"end": (-1, 0),
	"else": (-1, 1),
	"elif": (-1, 1)
}


def append_or_replace(in_list, in_index, in_item):
	if in_index < len(in_list):
		in_list[in_index] = in_item
	else:
		in_list.append(in_item)


def first_pass(in_file_lines):
	output_lines = []
	line_index = 0
	while line_index < len(in_file_lines):
		current_line = in_file_lines[line_index].strip()
		if re.search(r'^:\s{2,}', current_line) is not None:
			prev_line = in_file_lines[line_index - 1].strip()
			this_line = re.sub(r'^:\s{2,}', '', current_line).strip()
			in_file_lines[line_index - 1] = (prev_line + this_line).strip()
			in_file_lines[line_index] = ""
			line_index -= 1
			continue
		append_or_replace(output_lines, line_index, current_line)
		line_index += 1
	
	return output_lines


def clean_file(in_file_lines):
	global file_args
	f_lines = {fline.strip(): list() for fline in in_file_lines}
	i_lines = [fline.strip() for fline in in_file_lines]
	file_args = []

	n_lines = list(range(len(i_lines)))
	for tk, tt in line_replace.items():
		for ln in n_lines:
			ll = f"{i_lines[ln]}"
			if tt.test(ll):
				f_lines[ll].append(tk)

	o_lines = []
	c_il = 0
	for ll, tl in f_lines.items():
		rvl = ll
		n_in = 0
		for tk in tl:
			ci, ni = indent_change.get(tk, (0, 0))
			c_il += ci
			n_in += ni
			rvl = line_replace[tk].replace(rvl)
		rvl = ("\t" * c_il) + rvl
		o_lines.append(rvl)
		c_il += n_in

	o_lines = dict(enumerate(o_lines))
	while file_args:
		n_arg = file_args.pop().strip()
		for oln, ll in o_lines.items():
			for x in re.findall(f'[\\W]{n_arg}[\\W]', ll):
				n_x = x.replace(n_arg, n_arg.lower())
				ll = ll.replace(x, n_x)
			o_lines[oln] = ll

	return list(o_lines.values())


if __name__ == "__main__":
	converted_path = Path("converted")
	converted_path.mkdir(parents=True, exist_ok=True)
	
	parser = argparse.ArgumentParser(
		description='Format fortran in SLALib to nearly python'
	)
	parser.add_argument('in_file', nargs="+")
	parsed_args = parser.parse_args()
	
	with open("pypyslalib.py", "r") as fp:
		ff = fp.read()
		converted = set(re.findall("def (\\w+)", ff))
	
	inputted_files = set(parsed_args.in_file)
	file_list = set(x.stem.lower() for x in Path("pyslalib").glob("*.f"))
	inputted_files = inputted_files.union(file_list) - converted
	
	while inputted_files:
		in_fn = inputted_files.pop()
		if in_fn not in file_list:
			converted.add(in_fn)
			continue
		in_fpa = Path("pyslalib") / (in_fn + ".f")
		with open(in_fpa, "r") as if_fp:
			print("----------------------")
			print("converting ", in_fpa)
			input_file_lines = if_fp.read().splitlines()
			
			file_lines = first_pass(input_file_lines)
			file_lines = clean_file(file_lines)
			
			final_file_str = "\n".join(file_lines)
			final_file_str = autopep8.fix_code(
				final_file_str, options={'aggressive': 3}
			)
			
			ou_fpa = Path("converted") / (in_fn + ".py")
			with open(ou_fpa, "w") as cp_fp:
				cp_fp.write(final_file_str)
		converted.add(in_fn)
