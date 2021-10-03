import argparse
import glob
import os
import re
from pathlib import Path

line_replace = {
	"function": (
		r"SUBROUTINE (\w+) ?\(([ A-Za-z0-9,]+)\)",
		r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):"
	),
	"comment": (r"^\*+(.*)", r"# \g<1>"), "col_sign": (r"\s{2,}:\s{2,}", ""),
	"call": (r"CALL (\w+) ?\(([A-Za-z0-9,]+)\)", r"cls.\g<1>(\g<2>)"),
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
	"max": (r"MAX", r"np.maximum"), "end_file": ("END", ""), "implicit_none": ("IMPLICIT NONE", ""),
	"double_prec": (r"DOUBLE PRECISION (.+)", r""), "int": (r"INTEGER (.+)", r""),
	"double_constant": (r"(\d+)D(-?)(\d+)", r"\g<1>e\g<2>\g<3>")
}

if __name__ == "__main__":
	converted_path = Path("converted")
	converted_path.mkdir(parents=True, exist_ok=True)
	
	parser = argparse.ArgumentParser(description='Format fortran in SLALib to nearly python')
	parser.add_argument('in_file', nargs="+")
	args = parser.parse_args()
	
	for file in glob.iglob("pyslalib/*.f"):
		for file_arg in args.in_file:
			if file_arg.lower() in file.lower():
				with open(file, "r") as fp:
					file_lines = ""
					current_indent_level = 0
					for line in fp:
						line_indent_lv = current_indent_level
						new_line = line.strip()
						line_changed = False
						for val_type, val in line_replace.items():
							if re.search(val[0], new_line) is not None:
								line_changed = True
								new_line = re.sub(val[0], val[1], new_line)
								if val_type == "function":
									current_indent_level += 2
									match = re.search(r"sla_(\S+)\(", new_line)
									new_line = re.sub(r"sla_(\S+)\(", match.group(1).lower() + "(", new_line)
								if val_type in ("do", "if"):
									current_indent_level += 1
								if val_type == "end":
									current_indent_level -= 1
								if val_type in ("else", "elif"):
									line_indent_lv -= 1
								if val_type == "comment":
									new_line = re.sub(r"#\s+", r"# ", new_line)
									break
						file_lines += (line_indent_lv * "\t") + new_line + "\n"
					
					with open(os.path.join("converted", os.path.basename(file).split(".")[0] + ".py"), "w") as fp:
						fp.write(file_lines)
