import argparse
import json
import re
from boltons.setutils import IndexedSet
from pathlib import Path
from string import Template

import autopep8

file_args = IndexedSet()


def ind_set_extend(ind_set, in_iter):
	for item in in_iter:
		ind_set.add(item)


class Transformation:
	def __init__(self, f_test, py_replace):
		self.f_t = f_test
		self.py_r = py_replace
	
	def test(self, in_string):
		return re.search(self.f_t, in_string) is not None if in_string else False
	
	def replace(self, in_string):
		return re.sub(self.f_t, self.py_r, in_string)


class TemplateTransformation:
	def __init__(self, f_test, py_replace):
		self.f_t = Template(f_test)
		self.py_r = Template(py_replace)
		self.sub_dict = {}
	
	def __setitem__(self, key, value):
		self.sub_dict[key] = value
	
	def test(self, in_string):
		f_sub = self.f_t.substitute(**self.sub_dict)
		return re.search(f_sub, in_string) is not None if in_string else False
	
	def replace(self, in_string):
		f_sub = self.f_t.substitute(**self.sub_dict)
		r_sub = self.py_r.substitute(**self.sub_dict)
		return re.sub(f_sub, r_sub, in_string)


class DependT(Transformation):
	def replace(self, in_string):
		rv = super().replace(in_string)
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
		global pyf_info
		rv = super().replace(in_string)
		match = re.search(r"cls.(\w+)\((.+)\)", rv)
		if not match:
			return rv
		call_fn = match[1].lower()
		call_info = pyf_info[call_fn]
		call_in = call_info['in']
		call_out = call_info['out']
		
		call_fn_args = ",".join([x.strip() for x in match[2].split(",")[:len(call_in)]])
		rv = f"{','.join(call_out)} = cls.{call_fn}({call_fn_args})"
		
		new_dep = call_fn
		if new_dep not in converted:
			print("Adding Dependency ", new_dep)
			inputted_files.add(new_dep)
		return rv


class FuncT(Transformation):
	def replace(self, in_string):
		rv = super().replace(in_string)
		match = re.search(r"sla_(\w+)\((.+)\)", rv)
		fn_args = [x.strip() for x in match[2].split(",")]
		ind_set_extend(file_args, fn_args)
		fn_args = [x.lower() for x in fn_args]
		fn_args_str = ", ".join(fn_args)
		rv = re.sub(r"sla_(\w+)\((.+)\)", f'{match[1].lower()}({fn_args_str})', rv)
		return rv


fVar = r"[_A-Z][0-9A-Z_]+"

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
	"end": Transformation(r"END ?(DO|IF)?", r""),
	"elif": Transformation(r"ELSE IF", r"elif:"),
	"else": Transformation(r"ELSE", r"else:"),
	"if": Transformation(r"IF(.+)\sTHEN", r"if\g<1>:"),
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
	"arr_def": Transformation(r"DATA (\w+) / ([0-9,]+) /", "\g<1> = [\g<2>]"),
	"implicit_none": Transformation("IMPLICIT NONE", ""),
	"double_prec": Transformation(r"DOUBLE PRECISION (.+)", r""),
	"int": Transformation(r"INTEGER (.+)", r""),
	"double_constant": Transformation(
		r"(\d+)(?:D|E)(-?)\+?0*(\d+)", r"\g<1>e\g<2>\g<3>"
	),
	"function_call": CallT(r"sla_(\w+) ?\((.+)\)", r"cls.\g<1>(\g<2>)"),
	"file_return": TemplateTransformation(r"^sla_$cfn = (.+)", r"return \g<1>"),
	"sl_if": Transformation(r"IF\s?\((.+)\)\s?(.+)\s?=\s?(.+)", r"\g<2> = \g<3> if (\g<1>) else \g<2>"),
}

# (Local Change, Global Change)
indent_change = {
	"f_function": (0, 2),
	"function": (0, 2),
	"do": (0, 1),
	"if": (0, 1),
	"else": (-1, 0),
	"elif": (-1, 0),
	"end": (0, -1)
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


def clean_file(in_fn_args, fn_returns, in_file_name, in_file_lines):
	global file_args
	n_lines = len(in_file_lines)
	l_list = [f_line.strip() for f_line in in_file_lines]
	t_list = [list() for _ in range(n_lines)]
	file_args = IndexedSet(in_fn_args)
	
	line_replace['file_return']['cfn'] = in_file_name
	for ln in range(n_lines):
		ll = l_list[ln]
		for tk, tt in line_replace.items():
			if tt.test(ll):
				print("==================")
				print(ll, " passed testing for ", tk, " : ", tt.f_t)
				t_list[ln].append(tk)
	
	o_lines = []
	cur_ind_lvl = 0
	for l_num, l_trans in enumerate(t_list):
		tf_line = l_list[l_num]
		local_change = 0
		global_change = 0
		for tk in l_trans:
			l_cha, g_cha = indent_change.get(tk, (0, 0))
			local_change += l_cha
			global_change += g_cha
			tf_line = line_replace[tk].replace(tf_line)
		tf_line = ("\t" * (cur_ind_lvl + local_change)) + tf_line
		o_lines.append(tf_line)
		cur_ind_lvl += global_change
	
	r_str = ','.join(fn_returns)
	if 'sla_' not in r_str:
		o_lines.append(f"\t\treturn {','.join(fn_returns)}")
	
	o_lines = dict(enumerate(o_lines))
	non_ascii = r"[ *+=\-\)\(\\/]"
	while file_args:
		n_arg = file_args.pop().strip()
		for oln, ll in o_lines.items():
			for x in re.findall(f'{non_ascii}{n_arg}{non_ascii}', ll):
				n_x = x.replace(n_arg, n_arg.lower())
				ll = ll.replace(x, n_x)
			for x in re.findall(f'{non_ascii}{n_arg.upper()}{non_ascii}', ll):
				n_x = x.replace(n_arg.upper(), n_arg.lower())
				ll = ll.replace(x, n_x)
			o_lines[oln] = ll
	
	return list(o_lines.values())


def parse_pyf():
	if Path("parsed_pyf.json").exists():
		return json.load(open("parsed_pyf.json"))
	
	in_pyf = (Path("pyslalib") / "slalib.pyf").resolve().open("r")
	fn_name = ""
	fn_info = {}
	for line in in_pyf:
		if "function" in line or "subroutine" in line:
			if "end" not in line:
				fn_name = re.search(r"sla_(\w+)", line)[1]
				fn_info[fn_name] = {"in": [], "out": []}
		elif "::" in line:
			arg_def = False
			in_out_t = re.findall(r"intent ?\(([out,in ]+)\)", line)
			var_name = re.findall(r":: *(\w+)", line)
			if in_out_t:
				in_out_t = in_out_t[0]
			if var_name:
				var_name = var_name[0]
			if "out" in in_out_t:
				arg_def = True
				fn_info[fn_name]["out"].append(var_name)
			if "in" in in_out_t:
				arg_def = True
				fn_info[fn_name]["in"].append(var_name)
			if not arg_def:
				fn_info[fn_name]["in"].append(var_name)
	json.dump(fn_info, open("parsed_pyf.json", "w"), indent=2)
	return fn_info


if __name__ == "__main__":
	converted_path = Path("converted")
	converted_path.mkdir(parents=True, exist_ok=True)
	pyf_info = parse_pyf()
	parser = argparse.ArgumentParser(
		description='Format fortran in SLALib to nearly python'
	)
	parser.add_argument('in_file', nargs="+")
	parser.add_argument(
		'-l', '--limit', type=int, default=-1,
		help="Limit the number of files to convert, -1 is no limit"
	)
	
	parsed_args = parser.parse_args()
	
	with open("pypyslalib.py", "r") as fp:
		ff = fp.read()
		converted = IndexedSet(re.findall("def (\\w+)", ff))
	
	print(len(converted), " functions have been converted previously")
	
	inputted_files = IndexedSet(parsed_args.in_file)
	file_list = IndexedSet([x.stem.lower() for x in Path("pyslalib").glob("*.f")])
	inputted_files = inputted_files.union(file_list) - converted
	
	with open("internal_license.txt", "r") as fp:
		inline_license = [x.strip() for x in fp.read().strip().split(";;;;")]
	
	n_files_convert = 0
	while inputted_files:
		in_fn = inputted_files.pop()
		if in_fn not in file_list:
			converted.add(in_fn)
			continue
		in_fpa = Path("pyslalib") / f'{in_fn}.f'
		
		with open(in_fpa, "r") as if_fp:
			print("----------------------")
			print("converting ", in_fpa)
			file_lines = if_fp.read().splitlines()
			
			file_lines = first_pass(file_lines)
			file_lines = clean_file(pyf_info[in_fn]['in'], pyf_info[in_fn]['out'], in_fn, file_lines)
			
			final_file_str = "\n".join(file_lines)
			final_file_str = autopep8.fix_code(
				final_file_str, options={'aggressive': 3}
			)
			for il in inline_license:
				final_file_str = final_file_str.replace(il, "# ")
			ou_fpa = Path("converted") / f'{in_fn}.py'
			with open(ou_fpa, "w") as cp_fp:
				cp_fp.write(final_file_str)
		converted.add(in_fn)
		n_files_convert += 1
		if -1 < parsed_args.limit <= n_files_convert:
			break
