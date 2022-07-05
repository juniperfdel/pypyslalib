"""
PyPySLALib is licensed under GPLv3; 
PySLALib and SLALib are licensed under GPLv2;
The legal notices are below

---------------------------------------------------------------
PyPySLALib, a full python conversion of the astrometrics library SLALib
Copyright (C) 2022  Gregory Foote

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

---------------------------------------------------------------
f2py-generated wrappers for SLALIB
Copyright (C) 2010 Scott M. Ransom

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

---------------------------------------------------------------
SLALIB is a library of routines intended to make accurate and reliable positional-astronomy applications easier to write
Copyright (C) 1995 P.T.Wallace; Starlink; Rutherford Appleton Laboratory

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

---------------------------------------------------------------
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

import argparse
import json
import os
import re
import random
import sys

from pathlib import Path

import pyparsing as pp
from pyparsing import ParseException, pyparsing_common as ppc
from boltons.setutils import IndexedSet

all_fns_info = {}
converted = IndexedSet()
inputted_files = IndexedSet()
c_debug = False


class ConverterFailure(Exception):
    pass


class Singleton(object):
    _instance = None
    
    def __new__(cls, *args, **kwargs):
        if not isinstance(cls._instance, cls):
            cls._instance = object.__new__(cls, *args, **kwargs)
        return cls._instance


class FParser(Singleton):
    def __init__(self):
        self.var = ppc.identifier
        self.sci_num_sep = pp.oneOf("E D", caseless=True)
        self.num = pp.Combine(
            (pp.Word(f"+-{pp.nums}", pp.nums) + pp.Opt(pp.Literal(".") + pp.Opt(pp.Word(pp.nums)))) +
            pp.Opt(self.sci_num_sep + pp.Word(f"+-{pp.nums}", pp.nums))
        )
        self.var_num = (self.var | self.num)
        
        self.lpar = pp.Literal("(")
        self.rpar = pp.Literal(")")
        self.comma = pp.Literal(",")
        self.str_def = pp.oneOf("' \"")
        
        # PEMDAS
        self.pow = pp.oneOf("**")
        self.mul_op = pp.oneOf("* /")
        self.add_op = pp.oneOf("+ -")
        
        self.rel_op = pp.oneOf(".LT. .GT. .EQ. .LE. .GE. .NE.")
        
        self.comb_op = pp.oneOf(".AND. .OR. OR. AND.")
        
        self.equal_op = pp.Literal("=")
        
        self.expr = pp.Forward()
        self.atom = (
                (self.var + self.lpar + pp.OneOrMore(self.expr + pp.Opt(self.comma)) + self.rpar) |
                (self.num + pp.Opt(self.comma + self.expr)) |
                (self.var + pp.Opt(self.comma + self.expr)) |
                (self.lpar + pp.OneOrMore(self.expr + pp.Opt(self.comma)) + self.rpar) |
                pp.Combine(self.str_def + pp.SkipTo(self.str_def).leave_whitespace() + self.str_def)
        )
        
        self.factor = pp.Forward()
        self.factor << self.atom + pp.ZeroOrMore((self.pow + self.factor))
        self.term = pp.Opt(self.add_op) + self.factor + pp.ZeroOrMore((self.mul_op + self.factor))
        self.arith_expr = self.term + pp.ZeroOrMore((self.add_op + self.term))
        self.relational = self.arith_expr + pp.ZeroOrMore((self.rel_op + self.arith_expr))
        self.equal_terms = self.relational + pp.ZeroOrMore((self.comb_op + self.relational))
        self.expr << self.equal_terms + pp.ZeroOrMore((self.equal_op + self.equal_terms))
        
        CK = pp.CaselessKeyword
        self.subroutine = CK("SUBROUTINE") + self.var + pp.OneOrMore(self.expr + pp.Opt(","))
        self.d_function = (CK("REAL") | CK("DOUBLE PRECISION") | CK("INTEGER")) + \
                          CK("FUNCTION") + self.var + pp.OneOrMore(self.expr + pp.Opt(","))
        
        self.function = self.subroutine | self.d_function
        
        self.parameter = \
            (
                    CK("PARAMETER") +
                    pp.OneOrMore(self.expr + pp.Opt(","))
            )
        
        self.call = (
                CK("CALL") +
                self.var +
                self.lpar +
                pp.OneOrMore(self.expr + pp.Opt(",")) +
                self.rpar
        )
        
        self.i_for = (
                CK("DO") +
                self.expr
        )
        
        self.do_while = (
                CK("DO WHILE") +
                self.expr
        )
        
        self.end = \
            (CK("END IF") | CK("END DO"))
        
        self.el_if = (
            CK("ELSE IF") + pp.Opt(self.lpar) + self.expr + pp.Opt(self.rpar + self.expr) + pp.Opt(
            CK("THEN"))
        )
        
        self.else_t = CK("ELSE")
        
        self.if_t = (
                CK("IF") +
                pp.Opt(self.lpar) +
                self.expr +
                pp.Opt(self.rpar + self.expr) +
                pp.Opt(CK("THEN"))
        )
        
        self.arr = (
                CK("DATA") +
                pp.Opt(self.var_num) +
                pp.Literal("/") + pp.OneOrMore(self.expr + pp.Opt(",")) + pp.Literal("/")
        )
        
        self.int = (
                CK("INTEGER") +
                pp.OneOrMore(self.expr + pp.Opt(","))
        )
        
        self.char = (
                CK("CHARACTER") +
                self.var + pp.Literal(r"*(*)")
        )
        
        self.float = (
                (CK("REAL") | CK("DOUBLE PRECISION")) +
                pp.OneOrMore(self.expr + pp.Opt(","))
        )
        
        self.implicit = (
                CK("IMPLICIT") +
                pp.OneOrMore(self.expr + pp.Opt(","))
        )
        
        self.comment = (
                pp.Keyword(r"*") +
                pp.SkipTo(pp.LineEnd()).leave_whitespace() +
                pp.LineEnd()
        )
        
        self.stmt = (
                self.function.set_results_name("function") |
                self.call.set_results_name("call") |
                self.parameter.set_results_name("parameter") |
                self.do_while.set_results_name("while") |
                self.i_for.set_results_name("for") |
                self.end.set_results_name("end") |
                self.el_if.set_results_name("el_if") |
                self.else_t.set_results_name("else") |
                self.if_t.set_results_name("if") |
                self.arr.set_results_name("arr") |
                self.int.set_results_name("int") |
                self.float.set_results_name("float") |
                self.char.set_results_name("char") |
                self.implicit.set_results_name("implicit") |
                self.comment.set_results_name("comment") |
                self.expr.set_results_name("expr")
        )


class TObject:
    @staticmethod
    def test_token(in_token, in_str):
        try:
            in_token.parse_string(in_str)
            return True
        except ParseException:
            pass
        return False
    
    @staticmethod
    def get_arg_list(results, i_start):
        r_len = len(results)
        end_ind = 0
        
        lvl = 0
        commas = {0: [], 1: []}
        if results[i_start] != "(":
            raise ValueError("i_start must be the index of the starting parethesis for parsing!")
        for kk in range(i_start, r_len):
            if results[kk] == "(":
                lvl += 1
            if results[kk] == ")":
                lvl -= 1
                if lvl == 0:
                    end_ind = kk
                    break
            if results[kk] == ",":
                try:
                    commas[lvl].append(kk)
                except KeyError:
                    commas[lvl] = [kk]
        return end_ind, commas
    
    @staticmethod
    def add_dep(fn_name):
        global converted, inputted_files
        if fn_name not in converted:
            inputted_files.add(fn_name)
    
    @classmethod
    def test(cls, inst, i_type, results, ii):
        raise NotImplementedError
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        raise NotImplementedError


class TFn(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return (results[ii] in inst.fns) and (ii < inst.r_len - 1) and (results[ii+1] == "(")
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = inst.fns[results[ii]]


class TBadFn(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return (results[ii] in inst.bad_fns) and (ii < inst.r_len - 1) and (results[ii + 1] == "(")
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        clse_par, _ = cls.get_arg_list(results, ii+1)
        results[ii] = ""
        results[ii+1] = ""
        results[clse_par] = ""

class TKeyword(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return results[ii] in inst.kywds
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = inst.kywds[results[ii]]

class TNum(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return cls.test_token(inst.parse.num, results[ii])
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = results[ii].replace("D", "e").replace("E", "e")


class TForEq(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "for" and results[ii] == "="
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = " in range("
        for jj in range(ii + 1, inst.r_len, 2):
            if cls.test_token(inst.parse.num, results[jj]):
                results[jj] = str(int(results[jj]) - 1)
        results[inst.r_len - 1] += "):"


class TArrDef(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type in ("int", "float") and results[ii] == "(" and cls.test_token(inst.parse.var, results[ii - 1])
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        dtype_stmt = ", dtype=int" if i_type == "int" else ", dtype=float"

        inst.def_arrays.add(results[ii - 1])
        inst.arr_decl = True


        clse_par, cmmas_btwn = cls.get_arg_list(results, ii)
        if clse_par > 0:
            l_args = len(cmmas_btwn[1])
            if l_args > 1:
                results[ii] = " = np.zeros(("
                results[clse_par] = f"){dtype_stmt})"
            else:
                results[ii] = " = np.zeros("
                results[clse_par] = f"{dtype_stmt})"
            if (clse_par + 1) < inst.r_len and results[clse_par + 1] == ",":
                results[clse_par + 1] = ";"


class TVarDef(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type in ("int", "float") and ((ii == (inst.r_len - 1)) or (results[ii + 1] == ",")) and cls.test_token(
            inst.parse.var, results[ii]) and (not inst.arr_decl)
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        inst.def_vars.add(results[ii])
        if results[ii] not in inst.func_ins:
            results[ii] += " = 0" if i_type == "int" else " = 0."
            try:
                results[ii + 1] = ";"
            except IndexError:
                pass
        else:
            try:
                results[ii] = ""
                results[ii + 1] = ""
            except IndexError:
                pass


class TImpct(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "implicit"
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = ""


class TStrDef(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "char" and results[ii] == "*(*)"
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = "= \"\""

class TEnd(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "expr" and results[ii] == "END"
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = f"return {','.join(inst.func_outs)}"

class TFuncReturn(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "expr" and len(inst.func_name) > 0 and results[ii] == inst.func_name and results[ii + 1] == "="
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = "return"
        results[ii + 1] = ""

class TSingleLineIf(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "if" and results[ii] == ")" and (":" not in results[ii:])
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        inst.wc = 0
        rv = [*results[ii + 1:], *results[:ii + 1]]
        ei = 0
        for xi, xv in enumerate(rv):
            if xv == "=":
                ei = xi
                break
        rv[-1] += f" else {' '.join(rv[:ei])}"
        return rv, inst.r_len


class TSlaCheck(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "expr" and cls.test_token(inst.parse.var, results[ii]) and "sla_" in results[ii]
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        fn_name = results[ii].replace("sla_", "").lower()
        results[ii] = f"cls.{fn_name}"


class TArrUse(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return (results[ii] == "(") and cls.test_token(inst.parse.var, results[ii - 1]) and (
                results[ii - 1] in inst.def_arrays)
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        clse_par, cmmas_btwn = cls.get_arg_list(results, ii)
        results[ii] = "["
        if clse_par > 0:
            results[clse_par] = "]"


class TCallStmt(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return \
            i_type == "call" and \
            results[ii] == "(" and \
            cls.test_token(inst.parse.var, results[ii - 1]) and \
            results[ii - 1].startswith("sla_")
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        fn_name = results[ii - 1].replace("sla_", "").lower()
        cls.add_dep(fn_name)
        clse_par, cmmas_btwn = cls.get_arg_list(results, ii)
        tl_cs = cmmas_btwn[1]
        results[ii - 1] = f"cls.{fn_name}"
        
        in_args = len(all_fns_info[fn_name]['in'])
        t_args = len(tl_cs) + 1
        split_loc = clse_par if (t_args == 1) or (t_args <= in_args) else tl_cs[in_args - 1]
        
        out_args = results[split_loc+1:clse_par]
        call_and_in_args = results[:split_loc]
        after_call = results[clse_par:] if clse_par <= (inst.r_len - 1) else []
        call_and_in_args[ii - 1] = f" = {call_and_in_args[ii - 1]}"
        return out_args + call_and_in_args + after_call, inst.r_len


class TParamFix(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "parameter" and results[ii] == "("
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        clse_par, cmmas_btwn = cls.get_arg_list(results, ii)
        results[ii] = ""
        results[clse_par] = ""
        for ci in cmmas_btwn[1]:
            results[ci] = ";"
        return results, inst.r_len


class TUniNeg(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return \
            ii < (inst.r_len - 2) and \
            results[ii] == "-" and \
            (ii == 0 or (not cls.test_token(inst.parse.var_num, results[ii - 1]))) and \
            not cls.test_token(inst.parse.var_num, results[ii + 1])

    @classmethod
    def transform(cls, inst, i_type, results, ii):
        results[ii] = ""
        results[ii + 1] = f"-{results[ii + 1]}"


class TFuncDef(TObject):
    @classmethod
    def test(cls, inst, i_type, results, ii):
        return i_type == "function" and results[ii] == "(" and cls.test_token(inst.parse.var, results[ii - 1]) and \
               results[ii - 1].startswith("sla_")
    
    @classmethod
    def transform(cls, inst, i_type, results, ii):
        fn_name = results[ii - 1].replace("sla_", "").lower()
        cls.add_dep(fn_name)
        clse_par, cmmas_btwn = cls.get_arg_list(results, ii)
        in_args = len(all_fns_info[fn_name]['in'])
        if len(cmmas_btwn[1]) == 0:
            f_cl, arg_cl = clse_par, clse_par
        else:
            ind_a = in_args - 1
            f_cl, arg_cl = cmmas_btwn[1][0], cmmas_btwn[1][ind_a]
        
        inst.func_outs = ("".join(results[arg_cl + 1: clse_par])).split(",")
        inst.func_ins = ("".join(results[f_cl - 1:arg_cl + 1])).split(",")
        inst.func_ins = [x for x in inst.func_ins if len(x) > 0]
        inst.func_name = fn_name
        
        rv = [""]
        rv[0] = f"""
import numpy as np

class SLALib:
	@classmethod
	def {fn_name} (cls, {", ".join(inst.func_ins)}):"""
        return rv, inst.r_len


class FTransformer:
    stmts = {
        "FUNCTION": "",
        "SUBROUTINE": "",
        "PARAMETER": "",
        "CALL": "",
        "DO": "for ",
        "DO WHILE": "while ",
        "END IF": "",
        "END DO": "",
        "ELSE IF": "elif ",
        "ELSE": "else:",
        "IF": "if",
        "DATA": "",
        "INTEGER": "",
        "REAL": "",
        "CHARACTER": "",
        "DOUBLE PRECISION": "",
        "IMPLICIT": "",
        "*": "#"
    }
    
    kywds = {
        ".LT.": "<",
        ".GT.": ">",
        ".EQ.": "==",
        ".LE.": "<=",
        ".GE.": ">=",
        ".NE.": "!=",
        ".AND.": "and",
        ".OR.": "or",
        "AND.": "or",
        "OR.": "or",
        "THEN": ":",
        ".TRUE.": "True",
        ".FALSE.": "False",
        "REAL": "",
    }
    
    fns = {
        "ANINT": "np.rint",
        "SQRT": "np.sqrt",
        "SIGN": "np.sign",
        "NINT": "np.rint",
        "AINT": "np.trunc",
        "ABS": "np.abs",
        "DABS": "np.abs",
        "COS": "np.cos",
        "SIN": "np.sin",
        "TAN": "np.tan",
        "ATAN2": "np.arctan2",
        "MOD": "np.mod",
        "MAX": "np.maximum",
        "MIN": "np.minimum",
    }
    
    bad_fns = {
        "DBLE"
    }
    
    # (Local Change, Global Change)
    indent_change = {
        "function": (0, 2),
        "for": (0, 1),
        "while": (0, 1),
        "if": (0, 1),
        "else": (-1, 0),
        "el_if": (-1, 0),
        "end": (0, -1),
    }
    
    def __init__(self):
        self.def_arrays = set()
        self.def_vars = set()
        
        self.func_name = ""
        self.func_ins = set()
        self.func_outs = set()
        
        self.parse = FParser()
        
        self.ind_lvl = 0
        self.lc = 0
        self.wc = 0
        
        self.r_len = 0
        self.arr_decl = False
    
    def reset(self):
        self.__init__()
    
    def eval(self, i_string, parse_all=False):
        global c_debug
        if not i_string:
            return "\t" * self.ind_lvl
        
        try:
            results = self.parse.stmt.parse_string(i_string, parse_all=parse_all)
        except ParseException as e:
            raise ConverterFailure(f"Failed to convert {i_string}")
        
        i_type = list(results.as_dict().keys())[0]
        results = results.as_list()
        if results[0] in self.stmts:
            results[0] = self.stmts[results[0]]
        
        self.arr_decl = False
        
        self.r_len = len(results)
        self.lc, self.wc = self.indent_change.get(i_type, (0, 0))
        if c_debug:
            print("===================")
            print(i_string)
            print(i_type, " ;;;; ", results)
        for ii in range(self.r_len):
            
            # Transform Functions
            if TFn.test(self, i_type, results, ii):
                TFn.transform(self, i_type, results, ii)

            # Remove Converting/Fortran-specific Functions
            if TBadFn.test(self, i_type, results, ii):
                TBadFn.transform(self, i_type, results, ii)

            # Transform Keywords
            if TKeyword.test(self, i_type, results, ii):
                TKeyword.transform(self, i_type, results, ii)
            
            # Transform Numbers
            if TNum.test(self, i_type, results, ii):
                TNum.transform(self, i_type, results, ii)
            
            # Transform equal inside for into in range
            if TForEq.test(self, i_type, results, ii):
                TForEq.transform(self, i_type, results, ii)
            
            # Transform fortran array declarations
            if TArrDef.test(self, i_type, results, ii):
                TArrDef.transform(self, i_type, results, ii)
            
            # Recognize fortran variable declarations
            if TVarDef.test(self, i_type, results, ii):
                TVarDef.transform(self, i_type, results, ii)
            
            # Destroy IMPLICIT Statements
            if TImpct.test(self, i_type, results, ii):
                TImpct.transform(self, i_type, results, ii)
            
            # Recognize arrays being used in code
            if TArrUse.test(self, i_type, results, ii):
                TArrUse.transform(self, i_type, results, ii)
            
            # Build the return statement
            if TEnd.test(self, i_type, results, ii):
                TEnd.transform(self, i_type, results, ii)
            
            # Fix variables with sla_ in them
            if TSlaCheck.test(self, i_type, results, ii):
                TSlaCheck.transform(self, i_type, results, ii)
        
            # Fix the Fortran CHARACTER declaration
            if TStrDef.test(self, i_type, results, ii):
                TStrDef.transform(self, i_type, results, ii)
            
            # Fix Unary Negative sign
            if TUniNeg.test(self, i_type, results, ii):
                TUniNeg.transform(self, i_type, results, ii)
        
        # Second Pass to deal with overall restructuring of the code
        self.r_len = len(results)
        ii = 0
        while ii < self.r_len:
            # Transform CALL to classmethod calls
            if TCallStmt.test(self, i_type, results, ii):
                results, ii = TCallStmt.transform(self, i_type, results, ii)
            
            # Build the function definition
            elif TFuncDef.test(self, i_type, results, ii):
                results, ii = TFuncDef.transform(self, i_type, results, ii)
            
            # Fix single line if statements
            elif TSingleLineIf.test(self, i_type, results, ii):
                results, ii = TSingleLineIf.transform(self, i_type, results, ii)
            
            # Remove outermost parenthesis in parameter statements
            elif TParamFix.test(self, i_type, results, ii):
                results, ii = TParamFix.transform(self, i_type, results, ii)
            ii += 1
        
        results = ("\t" * (self.ind_lvl + self.lc)) + " ".join([xx.strip() for xx in results if len(xx) > 0])
        self.ind_lvl += self.wc

        if c_debug:
            print(results)
        return results


def remove_newline_op(in_file_lines):
    in_file_str = "\n".join(in_file_lines)
    in_file_str = re.sub(r"[\n\t\r ]+:[\n\t\r ]+", "", in_file_str).strip()
    return in_file_str.splitlines(False)


def clean_file(in_file_lines):
    f_trans = FTransformer()
    f_trans.reset()
    
    cleaned_file = []
    for ll in in_file_lines:
        if len(ll) < 1:
            cleaned_file.append("")
        rr = f_trans.eval(ll)
        cleaned_file.append(rr)
    return cleaned_file


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
    global all_fns_info, c_debug, converted, inputted_files
    
    converted_path = Path("converted")
    converted_path.mkdir(parents=True, exist_ok=True)
    all_fns_info = parse_pyf()
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
    
    inline_license = Path("internal_license.txt").read_text()
    
    n_files_convert = 0
    while inputted_files and parsed_args.limit:
        in_fn = inputted_files.pop()
        if (in_fn not in file_list) or (in_fn in converted):
            continue
        
        file_lines = (Path("pyslalib") / f"{in_fn}.f").resolve().read_text().splitlines()
        print("----------------------")
        print("converting ", in_fn)
        
        file_lines = remove_newline_op(file_lines)
        file_lines = clean_file(file_lines)
        
        final_file_str = "\n".join(file_lines)
        
        final_file_str = final_file_str.replace(inline_license, "\t\t#")
        final_file_str = final_file_str.replace(inline_license.replace("\t","    "), "        #")
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
