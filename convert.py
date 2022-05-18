import argparse
import json
import os
import re
import random
import sys
from pathlib import Path
from string import Template

import autopep8
from boltons.setutils import IndexedSet

file_args = IndexedSet()
known_arrs = IndexedSet()
converted = IndexedSet()
inputted_files = IndexedSet()
pyf_info = {}
c_debug = False


def ind_set_extend(ind_set, in_iter):
    for item in in_iter:
        ind_set.add(item)


class Transformation:
    def __init__(self, f_test, py_replace, priority_of=None):
        self.f_t = f_test
        self.py_r = py_replace
        self.priority_of = priority_of

    def __repr__(self):
        return f"<{type(self)}: ({self.f_t}, {self.py_r})>"

    def test(self, in_string):
        global c_debug
        in_string = in_string.strip()
        i_t = bool(re.search(self.f_t, in_string))
        if i_t and c_debug:
            print(self.f_t, " : ", in_string)
        return i_t

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
        in_string = in_string.strip()
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
        depends_ = [x for x in depends_ if (len(x) > 0) and (x not in converted)]
        for x in depends_:
            if x not in converted:
                inputted_files.add(x)
        return rv


class FuncDefT(Transformation):
    def replace(self, in_string):
        global pyf_info
        rv = super().replace(in_string)
        match = re.search(r"sla_(\w+)\((.+)\)", rv)
        fn_args = [x.strip() for x in match[2].split(",")]
        ind_set_extend(file_args, fn_args)

        call_fn = match[1].lower()
        call_info = pyf_info[call_fn]
        ind_set_extend(file_args, call_info["in"])
        ind_set_extend(file_args, call_info["out"])
        call_in = ["cls", *call_info["in"]]
        fn_args_str = ",".join([x.strip() for x in call_in])
        rv = re.sub(r"sla_(\w+)\((.+)\)", f"{call_fn}({fn_args_str})", rv)
        return rv


class FortranCallT(Transformation):
    def replace(self, in_string):
        global pyf_info
        rv = super().replace(in_string)
        match = re.search(r"cls.(\w+)\((.+)\)", rv)
        if not match:
            return rv
        call_fn = match[1].lower()
        if call_fn == "sleep":
            return f"time.sleep({match[2]})"
        call_fn = call_fn.replace("sla_", "")
        call_info = pyf_info[call_fn]
        call_in = call_info["in"]
        call_out = call_info["out"]

        call_fn_args = ",".join(
            [x.strip() for x in match[2].split(",")[: (len(call_in) + 1)]]
        )
        rv = f"{','.join(call_out)} = cls.{call_fn}({call_fn_args})"

        new_dep = call_fn
        if new_dep not in converted:
            print("Adding Dependency ", new_dep)
            inputted_files.add(new_dep)
        return rv


class FuncCallT(Transformation):
    def replace(self, in_string):
        rv = super().replace(in_string)
        match = re.search(r"cls\.(\w+)\(", rv)

        call_fn = match[1].lower()
        rv = re.sub(
            r"cls\.(\w+)\(",
            f"cls.{call_fn}(",
            rv,
        )
        return rv


class ArrDefT(Transformation):
    def replace(self, in_string):
        rv = super().replace(in_string)
        v_name, _ = rv.split("=", 1)
        v_name = v_name.strip()
        known_arrs.add(v_name)
        return rv


class ArrGetSetT(Transformation):
    def replace(self, in_string):
        rv = super().replace(in_string)
        if matches := re.findall(r"(\w+\s*)\[([\-\s\d]+)]", rv):
            for match in matches:
                arr_name = match[0].strip()
                arr_ind = int(match[1])
                rv = rv.replace(f"{match[0]}[{match[1]}]", f"{arr_name}[{arr_ind}]")
                known_arrs.add(arr_name)
        if "=" in rv:
            v_assign, v_value = rv.split("=", 1)
            v_value = re.sub(r"([ ,\b])\+(\d)", r"\g<1>\g<2>", v_value.strip())
            if "=" in v_value:
                v_value = re.sub(r",(?=\s+\w+\[[\-\s\d]+]\s+=)", "; ", v_value)

                t_l, n_l = v_value.split(";", 1)
                t_l = t_l.strip()
                v_value = f"{t_l};{n_l}"
            v_value = v_value.strip()
            v_assign = v_assign.strip()
            rv = f"{v_assign} = {v_value}"
        return rv


class LowerT(Transformation):
    def replace(self, in_string):
        m = re.search(self.f_t, in_string)
        rv = super().replace(in_string)
        for m_group in m.groups():
            rv = rv.replace(m_group, m_group.lower())
        return rv


class ParameterT(Transformation):
    def replace(self, in_string):
        rv = super().replace(in_string)
        if "(" in rv:
            rv = re.sub(r"\),?", "));", rv.replace("(", "= np.zeros(("))
        else:
            rv = rv.replace(",", ";")
        return rv


class BoolT(Transformation):
    def replace(self, in_string):
        return in_string.replace(".TRUE.", "True").replace(".FALSE.", "False")


line_replace = {
    "function": FuncDefT(
        r"SUBROUTINE ?(\w+) ?\(([ A-Za-z0-9,]+)\)",
        r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):",
    ),
    "f_function": FuncDefT(
        r".+FUNCTION (\w+) ?\(([ A-Za-z0-9,]+)\)",
        r"import numpy as np\n\n\nclass SLALib:\n\t@classmethod\n\tdef \g<1>(cls, \g<2>):",
    ),
    "dependencies": DependT(r"^\s*\*+\s*Called:(.+)", r"# Depends:\g<1>"),
    "comment": Transformation(r"^\*+\s*(.*)", r"# \g<1>"),
    "call": FortranCallT(r"CALL (\w+) ?\((.+)\)", r"cls.\g<1>(\g<2>)"),
    "do": Transformation(
        r"DO ([ \w]+)=([ \w\-]+),([ \w,\-]+)", r"for \g<1> in range(\g<2>, \g<3>):"
    ),
    "while": Transformation(r"DO WHILE ?(.+)", r"while \g<1>:"),
    "end": Transformation(r"END ?(DO|IF)?", r""),
    "elif": Transformation(r"ELSE IF(.+)(?:THEN)", r"elif \g<1>:", ("if", "else")),
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
    "anint": Transformation(r"(\W)ANINT(\W)", r"\g<1>np.rint\g<2>"),
    "dble": Transformation(r"(\W)DBLE(\W)", "\g<1>\g<2>"),
    "sqrt": Transformation(r"SQRT\(", r"np.sqrt("),
    "sign": Transformation(r"SIGN\(", r"np.sign("),
    "nint": Transformation(r"NINT\(", r"np.rint("),
    "abs": Transformation(r"ABS\(", r"np.abs("),
    "cos": Transformation(r"COS\(", r"np.cos("),
    "sin": Transformation(r"SIN\(", r"np.sin("),
    "tan": Transformation(r"TAN\(", r"np.tan("),
    "atan2": Transformation(r"ATAN2\(", r"np.arctan2("),
    "mod": Transformation(r"MOD\(", r"np.mod("),
    "max": Transformation(r"MAX\(", r"np.maximum("),
    "min": Transformation(r"MIN\(", "np.minimum("),
    "arr_set_2d": ArrGetSetT(
        r"(?<=DATA)(.*)\((\w+)\((\w+),([\-\s\d]+)\) ?, ?(\w+)=([\-\+\d,\s]+)\)\/([^\/]+)\/",
        r"\g<1> \g<2>[\g<4>] = [\g<7>]",
    ),
    "arr_def": ArrDefT(r"(?<=DATA) +(\w+) +/ ([^\/]+) /", r" \g<1> = [\g<2>]"),
    "arr_set_raw": ArrGetSetT(
        r"(?<=DATA) ?(\w+) ?\(([\-\s\d]+)\) ?\/ ?([^\/]+) ?\/",
        r" \g<1>[\g<2>] = \g<3>",
    ),
    "arr_set": ArrGetSetT(r"(\w+) ?\(([\-\s\d]+)\) ?= ?", r"\g<1>[\g<2>] = "),
    "arr_get": ArrGetSetT(
        r"(\w+) ?= ?(\w+) ?\(([A-Za-z0-9\- ]+)\)", r"\g<1> = \g<2>[\g<3>]"
    ),
    "implicit_none": Transformation("IMPLICIT NONE", ""),
    "parameter": ParameterT(
        r"(DOUBLE PRECISION|INTEGER|REAL|PARAMETER) ?\(?(.+)\)?\s*$",
        r"\g<2>",
        ("double_prec_bare", "int_bare", "real_bare", "type_redef"),
    ),
    "double_prec_bare": Transformation(r"DOUBLE PRECISION ([\w,]+)\s*$", r""),
    "int_bare": Transformation(r"INTEGER ([\w,]+)\s*$", r""),
    "real_bare": Transformation(r"REAL ([\w,]+)\s*$", r""),
    "data_def": Transformation(r"(DATA)", r""),
    "type_redef": LowerT(r"(INTEGER|REAL|FLOAT)\(", r"\g<1>("),
    "double_redef": Transformation("DBLE\(", "float("),
    "bool_constant": BoolT(r"\.(TRUE|FALSE)\.", r"\g<1>"),
    "double_constant": Transformation(
        r"(\d+)(?:D|E)(-?)\+?0*(\d+)", r"\g<1>e\g<2>\g<3>"
    ),
    "file_return": TemplateTransformation(r"^\s*sla_$cfn\s+=\s+(.+)", r"return \g<1>"),
    "function_call": FuncCallT(r"sla_(\w+) ?\((.+)\)", r"cls.\g<1>(\g<2>)"),
    "sl_if": Transformation(
        r"(?<!ELSE )IF\s?\((.+)\)\s?(.+)\s?=\s?(.+)",
        r"\g<2> = \g<3> if (\g<1>) else \g<2>",
    ),
}

# (Local Change, Global Change)
indent_change = {
    "f_function": (0, 2),
    "function": (0, 2),
    "do": (0, 1),
    "while": (0, 1),
    "if": (0, 1),
    "else": (-1, 0),
    "elif": (-1, 0),
    "end": (0, -1),
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
        if re.search(r"^:\s{2,}", current_line) is not None:
            prev_line = in_file_lines[line_index - 1].strip()
            this_line = re.sub(r"^:\s{2,}", "", current_line).strip()
            in_file_lines[line_index - 1] = (prev_line + this_line).strip()
            in_file_lines[line_index] = ""
            line_index -= 1
            continue
        append_or_replace(output_lines, line_index, current_line)
        line_index += 1

    return output_lines


def clean_fn_arr(in_olines):
    non_ascii = r"[ *+=\-\)\(\\/,_	]"
    find_args = r"\([_0-9A-Za-z]+\)"
    while known_arrs:
        n_arg = known_arrs.pop()
        for oln, ll in in_olines.items():
            search_str = f"{non_ascii}?{n_arg}{find_args}"
            for x in re.findall(search_str, ll, re.M):
                n_x = x.replace("(", "[").replace(")", "]")
                if matches := re.search(r"\[(\d+)\]", n_x):
                    arr_ind = int(matches[1]) - 1
                    n_x = n_x.replace(matches[0], f"[{arr_ind}]")
                ll = ll.replace(x, n_x)
            in_olines[oln] = ll


def clean_file(in_fn_args, fn_returns, in_file_name, in_file_lines):
    global file_args
    global known_arrs
    global c_debug
    n_lines = len(in_file_lines)
    l_list = [f_line.strip() for f_line in in_file_lines]
    t_list = [list() for _ in range(n_lines)]
    file_args = IndexedSet(in_fn_args)
    known_arrs = IndexedSet()

    line_replace["file_return"]["cfn"] = in_file_name.upper()
    for ln in range(n_lines):
        ll = l_list[ln]
        for tk, tt in line_replace.items():
            if tt.test(ll) and (
                tk != "function_call"
                or "f_function" not in t_list[ln]
                and "function" not in t_list[ln]
            ):
                t_list[ln].append(tk)

    if c_debug:
        print("====================================")

    o_lines = []
    cur_ind_lvl = 0
    for l_num, l_trans in enumerate(t_list):
        tf_line = l_list[l_num]
        destroy_ts = set()
        for tk in l_trans:
            if line_replace[tk].priority_of is not None:
                for p_of in line_replace[tk].priority_of:
                    destroy_ts.add(p_of)
        if c_debug:
            print(tf_line)
            print(l_trans)
            print(destroy_ts)
        local_change = 0
        global_change = 0
        for tk in l_trans:
            if tk in destroy_ts:
                continue
            l_cha, g_cha = indent_change.get(tk, (0, 0))
            local_change += l_cha
            global_change += g_cha
            tf_line = line_replace[tk].replace(tf_line)
        tf_line = ("\t" * (cur_ind_lvl + local_change)) + tf_line.strip()
        o_lines.append(tf_line)
        cur_ind_lvl += global_change

    r_str = ",".join(fn_returns)
    if "sla_" not in r_str:
        o_lines.append(f"\t\treturn {','.join(fn_returns)}")
        ind_set_extend(file_args, fn_returns)

    o_lines = dict(enumerate(o_lines))
    clean_fn_arr(o_lines)

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


def main():
    global file_args
    global known_arrs
    global pyf_info
    global converted
    global inputted_files
    global c_debug

    converted_path = Path("converted")
    converted_path.mkdir(parents=True, exist_ok=True)
    pyf_info = parse_pyf()
    parser = argparse.ArgumentParser(
        description="Format fortran in SLALib to nearly python"
    )
    parser.add_argument("in_file", nargs="*")
    parser.add_argument(
        "-l",
        "--limit",
        type=int,
        default=-1,
        help="Limit the number of files to convert, -1 is no limit",
    )
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--random", action="store_true")

    parsed_args = parser.parse_args()

    c_debug = parsed_args.debug

    with open("pypyslalib.py", "r") as fp:
        ff = fp.read()
        converted = IndexedSet(re.findall("def (\\w+)", ff))

    print(len(converted), " functions have been converted previously")
    inputted_files = IndexedSet(parsed_args.in_file)
    file_list = IndexedSet([x.stem.lower() for x in Path("pyslalib").glob("*.f")])
    if parsed_args.random:
        inputted_files.add(random.choice(file_list - converted))

    inputted_files = inputted_files - converted

    with open("internal_license.txt", "r") as fp:
        inline_license = [x.strip() for x in fp.read().strip().split(";;;;")]

    n_files_convert = 0
    while inputted_files and parsed_args.limit:
        in_fn = inputted_files.pop()
        if in_fn not in file_list:
            converted.add(in_fn)
            continue
        in_fpa = Path("pyslalib") / f"{in_fn}.f"

        with open(in_fpa, "r") as if_fp:
            print("----------------------")
            print("converting ", in_fpa)
            file_lines = if_fp.read().splitlines()

            file_lines = first_pass(file_lines)
            file_lines = clean_file(
                pyf_info[in_fn]["in"], pyf_info[in_fn]["out"], in_fn, file_lines
            )

            final_file_str = "\n".join(file_lines)
            final_file_str = autopep8.fix_code(
                final_file_str, options={"aggressive": 3}
            )
            for il in inline_license:
                final_file_str = final_file_str.replace(il, "# ")
            ou_fpa = Path("converted") / f"{in_fn}.py"
            with open(ou_fpa, "w") as cp_fp:
                cp_fp.write(final_file_str)
            try:
                os.system(f"{sys.executable} -m black {ou_fpa}")
            except OSError:
                pass
        converted.add(in_fn)
        n_files_convert += 1
        if -1 < parsed_args.limit <= n_files_convert:
            break


if __name__ == "__main__":
    main()
