import argparse
import glob
import os
import re
from pathlib import Path

import autopep8

line_replace = {
	"function": (
		r"SUBROUTINE (\w+) ?\(([ A-Za-z0-9,]+)\)",
		r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):"
	),
	"f_function": (
		r".+FUNCTION (\w+) ?\(([ A-Za-z0-9,]+)\)",
		r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):"
	),
	"dependencies": (r"^\s*\*+\s*Called:(.+)", r"# Depends:\g<1>"),
	"comment": (r"^\*+(.*)", r"# \g<1>"),
	"call": (r"CALL (\w+) ?\(([A-Za-z0-9,()\[\]._]+)\)", r"cls.\g<1>(\g<2>)"),
	"do": (r"DO (\S+)=(\S+),(\S+)", r"for \g<1> in (\g<2>, \g<3>):"),
	"end": (r"END (DO|IF)", r""), "eq_sign": ("=", "="),
	"elif": (r"ELSE IF", r"elif:"), "if": (r"IF(.+)\sTHEN", r"if\g<1>:"), "else": (r"ELSE", r"else:"),
	"sl_if": (r"IF \((.+)\)(.+)", r"if \g<1>: \g<2>"),
	"lt": (r" ?\.LT\. ?", r" < "), "gt": (r" ?\.GT\. ?", r" > "),
	"eq": (r" ?\.EQ\. ?", r" == "), "le": (r" ?\.LE\. ?", r" <= "),
	"ge": (r" ?\.GE\. ?", r" >= "), "neq": (r" ?\.NE\. ?", r" != "),
	"and": (r" ?\.AND\. ?", r" and "), "or": (r" ?\.OR\. ?", r" or "),
	"anint": (r"ANINT", r"np.rint"), "dble": (r"DBLE", ""),
	"sqrt": (r"SQRT", r"np.sqrt"), "sign": (r"SIGN", r"np.sign"),
	"nint": (r"NINT", r"np.rint"), "abs": (r"ABS", r"np.abs"),
	"cos": (r"COS\(", r"np.cos("), "sin": (r"SIN\(", r"np.sin("),
	"tan": (r"TAN\(", r"np.tan("), "atan2": (r"ATAN2\(", r"np.arctan2("),
	"mod": (r"MOD", r"np.mod"), "setter": (r"(\w+)\((\d+)\) = ", r"\g<1>[\g<2>] = "),
	"max": (r"MAX", r"np.maximum"), "min": (r"MIN", "np.minimum"), "end_file": ("END", ""),
	"implicit_none": ("IMPLICIT NONE", ""),
	"double_prec": (r"DOUBLE PRECISION (.+)", r""), "int": (r"INTEGER (.+)", r""),
	"double_constant": (r"(\d+)(?:D|E)(-?)\+?0*(\d+)", r"\g<1>e\g<2>\g<3>")
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


def clean_file(in_file_lines, in_line_replace, in_inputted_files):
	return_code = 1
	current_indent_level = 0
	output_file_lines = []
	b_s = r"\]\["
	for line in in_file_lines:
		line_indent_lv = current_indent_level
		new_line = line.strip()
		for val_type, val in in_line_replace.items():
			if val_type == "arguments":
				for arg_ in val:
					arg_regex = f"[ )(+-/*<>^{b_s}]*{arg_}[+-/*,= )(<>^{b_s}]+"
					match = re.search(arg_regex, new_line)
					if match is not None:
						new_line = re.sub(arg_regex, match.group(0).lower(), new_line)
			elif re.search(val[0], new_line) is not None:
				new_line = re.sub(val[0], val[1], new_line)
				if val_type in ("f_function", "function"):
					current_indent_level += 2
					match = re.search(r"sla_(\w+)\((.+)\)", new_line)
					fn_args = [x.strip() for x in match.group(2).split(",")]
					in_line_replace['arguments'].extend(fn_args.copy())
					fn_args = [x.lower() for x in fn_args]
					fn_args_str = ", ".join(fn_args)
					new_line = re.sub(r"sla_(\w+)\((.+)\)", match.group(1).lower() + "(" + fn_args_str + ")", new_line)
				if val_type == "call":
					match = re.search(r"sla_(\w+)\(", new_line)
					new_line = re.sub(r"sla_(\w+)\(", match.group(1).lower() + "(", new_line)
					print("Adding Dependency ", match.group(1).lower())
					in_inputted_files.append(match.group(1).lower())
					return_code = -1
				if val_type in ("do", "if"):
					current_indent_level += 1
				if val_type == "end":
					current_indent_level -= 1
				if val_type in ("else", "elif"):
					line_indent_lv -= 1
				if val_type == "comment":
					new_line = re.sub(r"#\s+", r"# ", new_line)
					break
				if val_type == "dependencies":
					depends_ = [x.replace("sla", "").replace("_", "").lower().strip() for x in new_line.split(":")[1].strip().split(",")]
					depends_ = [x for x in depends_ if len(x) > 0]
					
					in_inputted_files.extend(depends_)
					return_code = -1
		output_file_lines.append((line_indent_lv * "\t") + new_line + "\n")
	return return_code, output_file_lines


if __name__ == "__main__":
	converted_path = Path("converted")
	converted_path.mkdir(parents=True, exist_ok=True)
	
	parser = argparse.ArgumentParser(description='Format fortran in SLALib to nearly python')
	parser.add_argument('in_file', nargs="+")
	args = parser.parse_args()
	
	inputted_files = list(args.in_file)
	file_list = list(glob.iglob("pyslalib/*.f"))
	file_index = 0
	
	converted_files = set()
	while file_index < len(file_list):
		file = file_list[file_index]
		
		file_path_obj = Path(file)
		file_filename = file_path_obj.stem
		for file_arg in inputted_files:
			if (file not in converted_files) and (file_arg.lower() in file.lower()):
				this_line_replace = line_replace.copy()
				
				this_line_replace['file_return'] = (f"^sla_{file_filename.upper()} = (.+)", r"return \g<1>")
				this_line_replace['arguments'] = []
				with open(file, "r") as input_fortran_file:
					input_file_lines = input_fortran_file.read().splitlines()
				
				file_lines = first_pass(input_file_lines)
				rv_code, file_lines = clean_file(file_lines, this_line_replace, inputted_files)
				if rv_code < 0:
					file_index = -1
				
				final_file_str = "".join(file_lines)
				final_file_str = autopep8.fix_code(final_file_str, options={'aggressive': 3})
				with open(os.path.join("converted", os.path.basename(file).split(".")[0] + ".py"), "w") as converted_python_file_pointer:
					converted_python_file_pointer.write(final_file_str)
				
				inputted_files.remove(file_arg)
				converted_files.add(file)
		file_index += 1
