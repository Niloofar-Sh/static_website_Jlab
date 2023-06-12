import warnings
warnings.filterwarnings('ignore')
import os, base64,io
import sys
import xml.etree.ElementTree as ET
import rdflib
from rdflib import Graph, Namespace, URIRef
import re
from lxml import etree
import copy
import numpy as np
from numpyencoder import NumpyEncoder
import pandas as pd
import math
import operator as op
import xlsxwriter
import json
from SPARQLWrapper import SPARQLWrapper, JSON
import ipywidgets as widgets
from ipywidgets import interact, interact_manual, interactive, IntSlider, FloatSlider, Layout, FloatRangeSlider
from IPython.display import display, clear_output, Image, FileLink, HTML, Math, Javascript
from PIL import Image
from scipy.optimize import curve_fit, least_squares
import requests
from scipy.integrate import solve_ivp, LSODA
import matplotlib.pyplot as plt
from matplotlib import markers
import matplotlib.font_manager as font_manager
import matplotlib.colors
from pylab import rcParams
import matplotlib.image as mpimg
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import gzip
import csv
import glob
import urllib.request
import libcellml
from libcellml import Analyser, AnalyserModel, Component, Variable, Generator, GeneratorProfile, Importer, Model, Parser, Printer, Validator, Annotator, Units
import libsbml
from math import exp
import networkx as nx
from adjustText import adjust_text
import subprocess



def _dump_issues(source_method_name, logger):
    if logger.issueCount() > 0:
        print('The method "{}" found {} issues:'.format(source_method_name, logger.issueCount()))
        for i in range(0, logger.issueCount()):
            print('    - {}'.format(logger.issue(i).description()))


def parse_model(filename, strict_mode):
    cellml_file = open(filename)
    parser = Parser(strict_mode)
    model = parser.parseModel(cellml_file.read())
    # _dump_issues("parse_model", parser)
    return model


def print_model(model):
    printer = Printer()
    s = printer.printModel(model)
    return s


def validate_model(model):
    validator = Validator()
    validator.validateModel(model)
    # _dump_issues("validate_model", validator)
    return validator.issueCount()


def flatten_model(model, importer):
    flat_model = importer.flattenModel(model)
    return flat_model


def analyse_model(model):
    analyser = Analyser()
    a = analyser.analyseModel(model)
    # _dump_issues("analyse_model", analyser)
    return a


def resolve_imports(model, base_dir, strict_mode):
    importer = Importer(strict_mode)
    importer.resolveImports(model, base_dir)
    # _dump_issues("resolve_imports", importer)
    # if model.hasUnresolvedImports():
    #     print("unresolved imports?")
    # else:
    #     print("no unresolved imports.")
    return importer





def pValsFunc(singleSelection, floatSSget, dataset=False):
    # This section works for 'Template 3'
    if singleSelection.value == 'Template 3':

        solve = tuple()
        if dataset == False:
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ao$' + '  ($u$mol)'][0],)
        else:
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ao$' + '  ($u$mol)'][0],)

        bounds = [[], []]
        for i in range(6):
            bounds[0].append(0)
            bounds[1].append(np.inf)

        def func(X, p1, p2, p3, p4, p5, p6):
            Ai, Ao = X
            V_SS = (p1 * Ao - p2 * Ai) / (p3 + p4 * Ao + p5 * Ai * Ao + p6 * Ai)
            return V_SS

        if dataset == False:
            popt, pcov = curve_fit(func, solve,
                                   [el.value for el in floatSSget if el.description == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)
        else:
            popt, pcov = curve_fit(func, solve,
                                   [floatSSget[el] for el in floatSSget.columns if el == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)

        pVals = {}
        pVals[singleSelection.value] = {}
        for i, el in enumerate(popt):
            pVals[singleSelection.value]['p' + str(i + 1)] = el

        json_object = json.dumps(pVals, indent=4, sort_keys=True,
                                 separators=(', ', ': '), ensure_ascii=False, cls=NumpyEncoder)
        with open("./temporary files/pVals.json", "w") as outfile:
            outfile.write(json_object)

    # This section works for 'Template 4'
    if singleSelection.value == 'Template 4':
        solve = tuple()
        if dataset == False:
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bo$' + '  ($u$mol)'][0],)
        else:
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bo$' + '  ($u$mol)'][0],)

        bounds = [[], []]
        for i in range(10):
            bounds[0].append(0)
            bounds[1].append(np.inf)

        def func(X, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10):
            Ai, Ao, Bi, Bo = X
            V_SS = (p1 * Ao * Bi - p2 * Ai * Bo) / (
                        p3 * Ao + p4 * Ai + p5 * Bo + p6 * Bi + p7 * Ao * Bi + p8 * Ai * Ao + p9 * Ai * Bo + p10 * Bo * Bi)
            return V_SS

        if dataset == False:
            popt, pcov = curve_fit(func, solve,
                                   [el.value for el in floatSSget if el.description == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)
        else:
            popt, pcov = curve_fit(func, solve,
                                   [floatSSget[el] for el in floatSSget.columns if el == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)

        pVals = {}
        pVals[singleSelection.value] = {}
        for i, el in enumerate(popt):
            pVals[singleSelection.value]['p' + str(i + 1)] = el

        json_object = json.dumps(pVals, indent=4, sort_keys=True,
                                 separators=(', ', ': '), ensure_ascii=False, cls=NumpyEncoder)
        with open("./temporary files/pVals.json", "w") as outfile:
            outfile.write(json_object)

    # This section works for 'Template 7'
    if singleSelection.value == 'Template 7':
        solve = tuple()
        if dataset == False:
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bo$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ci$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Co$' + '  ($u$mol)'][0],)
        else:
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bo$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ci$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Co$' + '  ($u$mol)'][0],)

        bounds = [[], []]
        for i in range(13):
            bounds[0].append(0)
            bounds[1].append(np.inf)

        def func(X, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13):
            Ai, Ao, Bi, Bo, Ci, Co = X
            V_SS = (p1 * Ai * Bo * Co - p2 * Ao * Bi * Ci) / (
                        p3 * Ai + p4 * Ao + p5 * Bi + p6 * Bo + p7 * Ci + p8 * Co + p9 * Ai * Ao + p10 * Bi * Bo + p11 * Ci * Co + p12 * Ai * Bo * Co + p13 * Bi * Ci * Ao)
            return V_SS

        if dataset == False:
            popt, pcov = curve_fit(func, solve,
                                   [el.value for el in floatSSget if el.description == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)
        else:
            popt, pcov = curve_fit(func, solve,
                                   [floatSSget[el] for el in floatSSget.columns if el == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)

        pVals = {}
        pVals[singleSelection.value] = {}
        for i, el in enumerate(popt):
            pVals[singleSelection.value]['p' + str(i + 1)] = el

        json_object = json.dumps(pVals, indent=4, sort_keys=True,
                                 separators=(', ', ': '), ensure_ascii=False, cls=NumpyEncoder)
        with open("./temporary files/pVals.json", "w") as outfile:
            outfile.write(json_object)

    # This section works for 'Template 11'
    if singleSelection.value == 'Template 11':
        solve = tuple()
        if dataset == False:
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bo$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ci$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Co$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Di$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Do$' + '  ($u$mol)'][0],)
        else:
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bo$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ci$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Co$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Di$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Do$' + '  ($u$mol)'][0],)

        bounds = [[], []]
        for i in range(16):
            bounds[0].append(0)
            bounds[1].append(np.inf)

        def func(X, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16):
            Ai, Ao, Bi, Bo, Ci, Co, Di, Do = X
            V_SS = (p1 * Ai * Bo * Co * Do - p2 * Ao * Bi * Ci * Di) / (
                        p3 * Ai + p4 * Ao + p5 * Bi + p6 * Bo + p7 * Ci + p8 * Co + p9 * Di + p10 * Do + p11 * Ai * Ao + p12 * Bi * Bo + p13 * Ci * Co + p14 * Di * Do + p15 * Ai * Bo * Co * Do + p16 * Bi * Ci * Ao * Di)
            return V_SS

        if dataset == False:
            popt, pcov = curve_fit(func, solve,
                                   [el.value for el in floatSSget if el.description == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)
        else:
            popt, pcov = curve_fit(func, solve,
                                   [floatSSget[el] for el in floatSSget.columns if el == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)

        pVals = {}
        pVals[singleSelection.value] = {}
        for i, el in enumerate(popt):
            pVals[singleSelection.value]['p' + str(i + 1)] = el

        json_object = json.dumps(pVals, indent=4, sort_keys=True,
                                 separators=(', ', ': '), ensure_ascii=False, cls=NumpyEncoder)
        with open("./temporary files/pVals.json", "w") as outfile:
            outfile.write(json_object)

    # This section works for 'Template 5'
    if singleSelection.value == 'Template 5':

        solve = tuple()
        if dataset == False:
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == '$q_Bo$' + '  ($u$mol)'][0],)
            solve = solve + ([el.value for el in floatSSget if el.description == 'V_m' + '  (mV)'][0],)
        else:
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ai$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Ao$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bi$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == '$q_Bo$' + '  ($u$mol)'][0],)
            solve = solve + ([floatSSget[el] for el in floatSSget.columns if el == 'V_m' + '  (mV)'][0],)

        bounds = [[], []]
        for i in range(12):
            bounds[0].append(0)
            bounds[1].append(np.inf)

        def func(X, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12):
            Ai, Ao, Bi, Bo, V_m = X
            F = 96485
            R = 8.31
            T = 310
            if dataset == False:
                EXP = exp(F * V_m / (R * T))
            else:
                EXP = [exp(F * vm / (R * T)) for vm in V_m]

            V_SS = (p1 * Bo * Ao * EXP - p2 * Ai * Bi) / (
                        p11 * Bi * Ai * Ao * EXP + p12 * Bi * Bo * Ai * Ao * EXP + p3 * Bi * Bo + p4 * Bo * Ai + p5 * Bo * Ao * EXP + p6 * Bi * Ai + p7 * Bi * Ao * EXP + p8 * Bi * Bo * Ai + p9 * Bi * Bo * Ao * EXP + p10 * Bo * Ai * Ao * EXP);
            return V_SS

        if dataset == False:
            popt, pcov = curve_fit(func, solve,
                                   [el.value for el in floatSSget if el.description == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)
        else:
            popt, pcov = curve_fit(func, solve,
                                   [floatSSget[el] for el in floatSSget.columns if el == 'V_SS' + '  ($u$mol/s)'][0],
                                   maxfev=3000000, bounds=bounds)

        pVals = {}
        pVals[singleSelection.value] = {}
        for i, el in enumerate(popt):
            pVals[singleSelection.value]['p' + str(i + 1)] = el

        json_object = json.dumps(pVals, indent=4, sort_keys=True,
                                 separators=(', ', ': '), ensure_ascii=False, cls=NumpyEncoder)
        with open("./temporary files/pVals.json", "w") as outfile:
            outfile.write(json_object)
    return pVals





# Define a function to convert a infix expression to a MathML string (temporary solution with limitations 1. no support for units 2. some MathML elements defined in CellML2.0 are not supported)
def eq_to_mathml(infix, Vss, ode=False):
    if ode == False:
        preforumla = '<apply> <eq/> <ci>' + Vss + '</ci>'
    else:
        preforumla = '<apply> <eq/> <apply> <diff/> <bvar> <ci>' + 't' + '</ci> </bvar> <ci>' + 'X' + '</ci> </apply> \n'

    postformula = '</apply>'

    p = libsbml.parseL3Formula(infix)
    mathstr = libsbml.writeMathMLToString(p)
    # remove the <math> tags in the mathML string
    mathstr = mathstr.replace('<math xmlns="http://www.w3.org/1998/Math/MathML">', '')
    mathstr = mathstr.replace('</math>', '')
    mathstr = mathstr.replace('<?xml version="1.0" encoding="UTF-8"?>', '')
    # add left side of the equation
    mathstr1 = preforumla + mathstr + postformula
    mathstr2 = mathstr1.replace('\n', '')
    mathstr3 = mathstr2.replace('type="integer"', 'cellml:units="dimensionless"')
    mathstr = mathstr3.replace('<cn>', '<cn cellml:units="dimensionless">')

    return mathstr





def singleModelCellmlGenerator(modelUnits, savedpVals, button_modelName):
    # Setting up the units

    Concentration = Units("concentration_unit")
    if modelUnits['Concentration'] == 1:
        Concentration.addUnit("mole")  # reference unit and built-in prefix
    else:
        Concentration.addUnit("mole", modelUnits['Concentration'])

    flowRate = Units("flow_rate_unit")
    flowRate.addUnit("concentration_unit", 1)  # reference unit, prefix, exponent, multiplier
    if modelUnits['Flow rate'] == 1:
        flowRate.addUnit('second', 1, -1)
    else:
        flowRate.addUnit('second', modelUnits['Flow rate'], -1)

    membraneVoltageUnit = Units("mV")
    membraneVoltageUnit.addUnit('volt', 'milli')

    C_per_mol = Units("C_per_mol")
    C_per_mol.addUnit("coulomb", 1)  # reference unit, prefix, exponent, multiplier
    C_per_mol.addUnit('mole', 1, -1)

    J_per_K_per_mol = Units("J_per_K_per_mol")
    J_per_K_per_mol.addUnit("joule", 1)  # reference unit, prefix, exponent, multiplier
    J_per_K_per_mol.addUnit("kelvin", 1, -1)
    J_per_K_per_mol.addUnit("mole", 1, -1)

    # Setting up the maths

    if singleSelection.value == 'Template 3':
        equation = '(p1* Ao - p2* Ai)/ (p3 + p4*Ao + p5*Ai*Ao + p6*Ai)'
    if singleSelection.value == 'Template 4':
        equation = '(p1*Ao*Bi - p2*Ai*Bo)/ (p3*Ao + p4*Ai + p5*Bo + p6*Bi + p7*Ao*Bi + p8*Ai*Ao + p9*Ai*Bo + p10*Bo*Bi)'
    if singleSelection.value == 'Template 7':
        equation = '(p1*Ai*Bo*Co - p2*Ao*Bi*Ci)/(p3*Ai + p4*Ao + p5*Bi + p6*Bo + p7*Ci + p8*Co + p9*Ai*Ao + p10*Bi*Bo + p11*Ci*Co + p12*Ai*Bo*Co + p13*Bi*Ci*Ao )'
    if singleSelection.value == 'Template 11':
        equation = '(p1*Ai*Bo*Co*Do - p2*Ao*Bi*Ci*Di)/(p3*Ai + p4*Ao + p5*Bi + p6*Bo + p7*Ci + p8*Co + p9*Di + p10*Do + p11*Ai*Ao + p12*Bi*Bo + p13*Ci*Co + p14*Di*Do + p15*Ai*Bo*Co*Do + p16*Bi*Ci*Ao*Di )'
    if singleSelection.value == 'Template 5':
        equation = '(p1*Bo*Ao*exp(F*V_m/(R*T))-p2*Ai*Bi)/(p11*Bi*Ai*Ao*exp(F*V_m/(R*T))+p12*Bi*Bo*Ai*Ao*exp(F*V_m/(R*T))+p3*Bi*Bo+p4*Bo*Ai+p5*Bo*Ao*exp(F*V_m/(R*T))+p6*Bi*Ai+p7*Bi*Ao*exp(F*V_m/(R*T))+p8*Bi*Bo*Ai+p9*Bi*Bo*Ao*exp(F*V_m/(R*T))+p10*Bo*Ai*Ao*exp(F*V_m/(R*T)))'

    # for variable in valueRequired[singleSelection.value]:
    #     if variable not in [x for x in valueRequired1[singleSelection.value]]:
    mathBody = eq_to_mathml(equation, 'V_SS', ode=False)

    fakeODE = eq_to_mathml('0', 'X', ode=True)

    # Generate the CellML model, components, variables,...

    model = libcellml.Model()
    model.setName(button_modelName)

    model.addUnits(Concentration)
    model.addUnits(flowRate)
    model.addUnits(membraneVoltageUnit)
    model.addUnits(C_per_mol)
    model.addUnits(J_per_K_per_mol)

    # Create a new component
    component = libcellml.Component()
    component.setName("main")

    # Add the component to the model
    model.addComponent(component)

    v = libcellml.Variable()
    v.setName('t')
    v.setUnits('second')
    v.setInterfaceType("public")
    component.addVariable(v)

    v = libcellml.Variable()
    v.setName('T')
    v.setUnits('kelvin')
    v.setInitialValue(310)
    v.setInterfaceType("public")
    component.addVariable(v)

    v = libcellml.Variable()
    v.setName('F')
    v.setUnits('C_per_mol')
    v.setInitialValue(96485)
    v.setInterfaceType("public")
    component.addVariable(v)

    v = libcellml.Variable()
    v.setName('R')
    v.setUnits('J_per_K_per_mol')
    v.setInitialValue(8.31)
    v.setInterfaceType("public")
    component.addVariable(v)

    # Create variables
    # Concentration variables here,
    for variable in valueRequired[singleSelection.value]:
        if variable in [x for x in valueRequired1[singleSelection.value]]:
            v = Variable()
            v.setName(variable)
            v.setUnits(Concentration)
            v.setInitialValue(
                [x.value for x in floatSSget if x.description == '$q_{}$'.format(variable) + '  ($u$mol)'][0])
            v.setInterfaceType("public")
            component.addVariable(v)
        else:
            if variable == 'V_SS':
                v = Variable()
                v.setName(variable)
                v.setUnits(flowRate)
                v.setInterfaceType("public")
                component.addVariable(v)
            if variable == 'V_m':
                v = Variable()
                v.setName(variable)
                v.setUnits(membraneVoltageUnit)
                v.setInitialValue([x.value for x in floatSSget if x.description == variable + '  (mV)'][0])
                v.setInterfaceType("public")
                component.addVariable(v)

    # f = open("./savedModels/savedpVals.json")
    # savedpVals = json.load(f)
    for variable in savedpVals[button_modelName]:
        v = Variable()
        v.setName(variable)
        v.setUnits('dimensionless')
        v.setInitialValue(savedpVals[button_modelName][variable])
        v.setInterfaceType("public")
        component.addVariable(v)

    # Setting up a fake ODE to run the code
    v = Variable()
    v.setName('X')
    v.setUnits('dimensionless')
    v.setInitialValue('0')
    v.setInterfaceType("public")
    component.addVariable(v)

    # Assign IDs
    annotator = Annotator()
    annotator.setModel(model)

    annotator.clearAllIds()
    annotator.assignAllIds()

    # Create the equation
    math_header = '<math xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/2.0#">'
    math_footer = '</math>'

    component.setMath(math_header)
    component.appendMath(mathBody)
    component.appendMath(fakeODE)
    component.appendMath(math_footer)

    printer = libcellml.Printer()

    writeCellML(model, printer, button_modelName)

    return model





# Write a model to cellml file, input: directory, model, output: cellml file
def writeCellML(model, printer, button_modelName):
    #     Write the serialised string to a file.
    write_file = open(button_modelName + ".cellml", "w")
    write_file.write(printer.printModel(model))
    print("The model has been printed to: {}.cellml".format(model.name()))
    return





def xmlAnnot(model, savedAnnotations):
    # Define namespaces
    rdf = Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    bqbiol = Namespace("http://biomodels.net/biology-qualifiers/")

    # Create RDF graph
    g = Graph()
    g.bind("rdf", rdf)
    g.bind("bqbiol", bqbiol)

    for compNum in range(model.componentCount()):
        for varNum in range(model.component(compNum).variableCount()):
            if model.component(compNum).variable(varNum).name() in ['Ao', 'Ai', 'Bo', 'Bi', 'Co', 'Ci', 'Do', 'Di']:
                to_be_annotated_var = model.component(compNum).variable(varNum).id()
                # Define RDF resources for the Concentrations
                Q1 = URIRef(button_modelName.value.strip('\n') + ".cellml" + '#' + to_be_annotated_var)
                opb = URIRef("https://identifiers.org/opb/" + [x.replace('isVersionOf$', '') for x in
                                                               savedAnnotations[button_modelName.value.strip('\n')][
                                                                   model.component(compNum).variable(varNum).name()] if
                                                               'isVersionOf$' in x][0])
                for x in savedAnnotations[button_modelName.value.strip('\n')][
                    model.component(compNum).variable(varNum).name()]:
                    if 'CHEBI' in x:
                        entity = URIRef("https://identifiers.org/chebi/" + [x.replace('entity$', '') for x in
                                                                            savedAnnotations[
                                                                                button_modelName.value.strip('\n')][
                                                                                model.component(compNum).variable(
                                                                                    varNum).name()] if 'entity$' in x][
                            0])
                    if 'GO' in x:
                        entity = URIRef("http://purl.obolibrary.org/obo/" + [x.replace('entity$', '') for x in
                                                                             savedAnnotations[
                                                                                 button_modelName.value.strip('\n')][
                                                                                 model.component(compNum).variable(
                                                                                     varNum).name()] if 'entity$' in x][
                            0])
                    if 'hasSourceParticipant$' in x:
                        sourceSink = URIRef("http://bime.uv.edu/semsim/" +
                                            [x.replace('hasSourceParticipant$', '') for x in
                                             savedAnnotations[button_modelName.value.strip('\n')][
                                                 model.component(compNum).variable(varNum).name()] if
                                             'hasSourceParticipant$' in x][0])
                    if 'hasSinkParticipant$' in x:
                        sourceSink = URIRef("http://bime.uv.edu/semsim/" +
                                            [x.replace('hasSinkParticipant$', '') for x in
                                             savedAnnotations[button_modelName.value.strip('\n')][
                                                 model.component(compNum).variable(varNum).name()] if
                                             'hasSinkParticipant$' in x][0])

                fma = URIRef("http://identifiers.org/fma/" + [x.replace('isPartOf$', '') for x in
                                                              savedAnnotations[button_modelName.value.strip('\n')][
                                                                  model.component(compNum).variable(varNum).name()] if
                                                              'isPartOf$' in x][0])

                # Add RDF triples to graph
                g.add((Q1, bqbiol.isVersionOf, opb))
                g.add((Q1, bqbiol.isVersionOf, entity))
                g.add((Q1, bqbiol.isPartOf, fma))
                g.add((Q1, bqbiol.sourceSink, sourceSink))

            if model.component(compNum).variable(varNum).name() == 'V_SS':
                to_be_annotated_var = model.component(compNum).variable(varNum).id()
                # Define RDF resources for the Concentrations
                Q1 = URIRef(button_modelName.value.strip('\n') + ".cellml" + '#' + to_be_annotated_var)
                opb = URIRef("https://identifiers.org/opb/" + [x.replace('isVersionOf$', '') for x in
                                                               savedAnnotations[button_modelName.value.strip('\n')][
                                                                   model.component(compNum).variable(varNum).name()] if
                                                               'isVersionOf$' in x][0])

                g.add((Q1, bqbiol.isVersionOf, opb))

    # Serialize RDF graph as XML file
    xml_string = g.serialize(format='xml')

    # and save to a file
    with open('./{}.xml'.format(button_modelName.value.strip('\n')), 'w') as f:
        f.write(g.serialize(format='xml'))





def convertToCellml2(cellml_file, mainDir, newCellmlName):
    cellml_strict_mode = False
    if len(sys.argv) > 2:
        strict_mode = sys.argv[2]
        if strict_mode == 'false':
            cellml_strict_mode = False

    # if cellml_strict_mode:
    #     print('  Parsing files in STRICT mode (only CellML 2.0 models accepted)')
    # else:
    #     print('  Parsing files in NON-STRICT mode (any CellML models accepted)')

    model = parse_model(cellml_file, cellml_strict_mode)
    if validate_model(model) > 0:
        exit(-1)

    importer = resolve_imports(model, mainDir, cellml_strict_mode)
    if model.hasUnresolvedImports():
        print("unresolved imports?")
        exit(-2)

    if validate_model(model) > 0:
        print('Validation issues found')
        exit(-3)

    # print('Model was parsed, resolved, and validated without any issues.')

    # need a flattened model for analysing
    flat_model = flatten_model(model, importer)
    if validate_model(model) > 0:
        # print('Validation issues found in flattened model')
        exit(-4)

    # print('Model was flattened without any issues.')

    # this will report any issues that come up in analysing the model to prepare for code generation
    analysed_model = analyse_model(flat_model)

    annotator = Annotator()
    annotator.setModel(model)

    annotator.clearAllIds()
    annotator.assignAllIds()
    model_string = print_model(model)

    # and save the updated model to a new file - note, we need the model filename for making our annotations later
    model_file = './savedModels/PMR/{}.cellml'.format(newCellmlName)
    with open(model_file, 'w') as f:
        f.write(model_string)
    return model





def pmrSearching(name):
    sparqlendpoint = 'https://models.physiomeproject.org/pmr2_virtuoso_search'
    sparql = SPARQLWrapper(sparqlendpoint)

    def search_entity(terms):
        query = """SELECT ?graph ?Model_entity WHERE {{ GRAPH ?graph {{ ?Model_entity ?p ?o FILTER REGEX(LCASE(STR(?Model_entity)), '{terms}')}}}}""".format(
            terms=terms)
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        graphs = sparql.query().convert()
        return graphs

    def search_model(terms):
        terms = terms.lower()
        entities = search_entity(terms)
        model = set()
        for entity in entities['results']['bindings']:
            workspace = entity['graph']['value']
            cellml = entity['Model_entity']['value'].split('#')[0]
            if not cellml.startswith('http') and terms in cellml.lower():
                model.update([workspace + '/rawfile/HEAD/' + cellml])
        return list(model)

    pmrModel = search_model(name)

    return pmrModel





def RDFpmrSearching(name):
    sparqlendpoint = 'https://models.physiomeproject.org/pmr2_virtuoso_search'
    sparql = SPARQLWrapper(sparqlendpoint)

    def search_entity(terms):
        query = """SELECT ?graph ?s ?p ?o WHERE {{ GRAPH ?graph {{ ?s ?p ?o FILTER REGEX(LCASE(STR(?s)), '{terms}')}}}}""".format(
            terms=terms)
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        graphs = sparql.query().convert()
        return graphs

    def search_model(terms):
        terms = terms.lower()
        entities = search_entity(terms)
        results = {}
        for entity in entities['results']['bindings']:
            workspace = entity['graph']['value']
            if workspace not in results:
                results[workspace] = []
            results[workspace] += [(entity['s']['value'], entity['p']['value'], entity['o']['value'])]
        return results

    return search_model(name)





def q_to_mathml(infix, speciesNum, ode=False):
    if ode == False:
        preforumla = '<apply> <eq/> <ci>' + 'q' + str(speciesNum) + '</ci>'
    else:
        preforumla = '<apply> <eq/> <apply> <diff/> <bvar> <ci>' + 't' + '</ci> </bvar> <ci>' + 'q' + str(
            speciesNum) + '</ci> </apply> \n'

    postformula = '</apply>'
    p = libsbml.parseL3Formula(infix)
    mathstr = libsbml.writeMathMLToString(p)
    # remove the <math> tags in the mathML string
    mathstr = mathstr.replace('<math xmlns="http://www.w3.org/1998/Math/MathML">', '')
    mathstr = mathstr.replace('</math>', '')
    mathstr = mathstr.replace('<?xml version="1.0" encoding="UTF-8"?>', '')
    # add left side of the equation
    mathstr1 = preforumla + mathstr + postformula
    mathstr2 = mathstr1.replace('\n', '')
    mathstr3 = mathstr2.replace('type="integer"', 'cellml:units="dimensionless"')
    mathstr = mathstr3.replace('<cn>', '<cn cellml:units="dimensionless">')

    return mathstr





def pythonCreator(mainDir, new_cellml_file):
    # if __name__ == '__main__':

    # STEP 1
    # Parse the model from a CellML file.

    # Create a libCellML Parser, and use it to parse the fileContents
    # string and convert it into a CellML Model structure.
    read_file = open('{}{}.cellml'.format(mainDir, new_cellml_file), "r")
    parser = Parser()
    model = parser.parseModel(read_file.read())

    # STEP 2
    # Resolve any import dependencies (if present) in the model.

    if (model.hasUnresolvedImports()):
        # Create an Importer instance.
        importer = Importer()

        # Submit the model to the importer and the absolute location
        # against which the import reference paths will be resolved.
        importer.resolveModelImports(model, "resources/")
        print_issues_to_terminal(importer)

        # Print a list of sources that this model requires. This list will
        # be empty after the model has been flattened.
        print_import_dependencies(model)

        # Retrieve a "flattened" (ie: import-free) model from the importer,
        # and use it to over-write the current model.
        model = importer.flattenModel(model)

    # STEP 3
    # Validate the model: check for syntactic and semantic errors.

    # Create a Validator instance and pass the model for checking.
    validator = Validator()
    validator.validateModel(model)

    # STEP 4
    # Analyse a model: check for mathematical and modelling errors.
    analyser = Analyser()
    analyser.analyseModel(model)

    if validate_model(model) > 0:
        print('Validation issues found')
        exit(-1)

    if model.hasUnresolvedImports():
        print("unresolved imports?")
        exit(-2)

    print('Model was parsed, resolved, and validated without any issues.')

    # STEP 5
    # Generate runnable code in other language formats for this model.

    # Create a Generator instance.  Note that by default this is the C language.
    generator = Generator()
    generator.setModel(analyser.model())
    generator.interfaceCode()
    generator.implementationCode()

    # Pass the generator the analysed model for processing.
    # generator.processModel(analyser.model())

    # STEP 6

    # If required, change the generator profile to Python and reprocess the model.
    profile = GeneratorProfile(GeneratorProfile.Profile.PYTHON)
    generator.setProfile(profile)
    # generator.processModel(model)

    # Retrieve and write the implementation code (*.py) to a file.
    write_file = open('{}.py'.format(new_cellml_file), "w")
    write_file.write(generator.implementationCode())
    write_file.close()

    # END





def PMRmodelComposition():
    mainDir = os.path.dirname("./savedModels/PMR/")
    folder_path = "./savedModels/PMR"
    extensions = (".cellml")

    global newCellmlNames
    cellml_files = []
    newCellmlNames = []

    # Loop through files in folder
    for file in os.listdir(folder_path):
        # Check if the file has a specific extension
        if file.endswith(extensions):
            # Add the file to the list
            cellml_files.append(os.path.join(folder_path, file))

    for file in cellml_files:
        filename = os.path.basename(file).split('.')[0]
        newCellmlNames.append(filename)

    newModels = {}
    for cellml_file, newCellmlName in zip(cellml_files, newCellmlNames):
        newModels[newCellmlName] = convertToCellml2(cellml_file, mainDir, newCellmlName)

    # global savedAnnotations
    savedAnnotations = {}

    for key in xmlGrouped:
        newKey = key.replace('.cellml', '')
        savedAnnotations[newKey] = xmlGrouped[key]

    global Species
    Species = []
    for transporterName in savedAnnotations:
        for entity in savedAnnotations[transporterName]:
            for predicate, obj in savedAnnotations[transporterName][entity]:
                if 'https://identifiers.org/opb/OPB_00425' in obj:  # identifying the concentration of chemicals
                    Species.append([transporterName, savedAnnotations[transporterName][entity]])

    global speciesNoDuplicate
    speciesNoDuplicate = []

    for i, element in enumerate(Species):
        # remove 'sink' and 'source' strings from the element and sort the list
        clean_element = sorted([x for x in element[1] if
                                x[1] not in ["http://bime.uv.edu/semsim/Source", "http://bime.uv.edu/semsim/Sink"]])
        # check if the element is unique
        is_unique = True
        for j in range(i + 1, len(Species)):
            if sorted([x for x in Species[j][1] if x[1] not in ["http://bime.uv.edu/semsim/Source",
                                                                "http://bime.uv.edu/semsim/Sink"]]) == clean_element:
                is_unique = False
                break
        if is_unique:
            speciesNoDuplicate.append(clean_element)

    for i, species in enumerate(speciesNoDuplicate):
        speciesNoDuplicate[i] = speciesNoDuplicate[i] + ['q' + str(i)]

    V = []
    for transporterName in savedAnnotations:
        for entity in savedAnnotations[transporterName]:
            for predicate, obj in savedAnnotations[transporterName][entity]:
                if 'https://identifiers.org/opb/OPB_00592' in obj:  # identifying the SS fluxes
                    V.append([transporterName, savedAnnotations[transporterName][entity], 'V_SS'])

    M = np.zeros((len(speciesNoDuplicate), len(V)))
    coefficient = 0
    for i, species in enumerate(speciesNoDuplicate):
        for element in Species:
            if all(el in element[1] for el in species[:-1]):
                for y in element[1]:
                    if 'http://bime.uv.edu/semsim/Sink' in y[1]:
                        coefficient = +1
                    if 'http://bime.uv.edu/semsim/Source' in y[1]:
                        coefficient = -1
                    for j, fluxinfo in enumerate(V):
                        if element[0] == fluxinfo[0]:
                            M[i][j] = coefficient

    odeRightSide = {}
    for i in range(M.shape[0]):
        odeRightSide[i] = ''
        for j in range(len(V)):
            if M[i][j] > 0:
                odeRightSide[i] += '+' + str(M[i][j]) + '*' + V[j][0].replace('.cellml', '') + '_' + V[j][2]
            if M[i][j] < 0:
                odeRightSide[i] += str(M[i][j]) + '*' + V[j][0].replace('.cellml', '') + '_' + V[j][2]

    cellMLRef = {}
    for modelName in newModels:
        for compNum in range(newModels[modelName].componentCount()):
            for varNum in range(newModels[modelName].component(compNum).variableCount()):
                for flux in V:
                    if flux[0] == modelName:
                        for ID, value in savedAnnotations[modelName].items():
                            if all(el in flux[1][0] for el in value[0]) and ID == newModels[modelName].component(
                                    compNum).variable(varNum).id():
                                cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(
                                    varNum).name()] = []
                                cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(
                                    varNum).name()].append(
                                    newModels[modelName].component(compNum).variable(varNum).units())
                                cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(
                                    varNum).name()].append(
                                    newModels[modelName].component(compNum).variable(varNum).initialValue())
                                cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(
                                    varNum).name()].append(
                                    newModels[modelName].component(compNum).variable(varNum).id())

                for speciesNum, species in enumerate(speciesNoDuplicate):
                    for value in savedAnnotations[modelName].values():
                        if all(el in value[:-1] for el in species[:-1]) and value[-1] == newModels[modelName].component(
                                compNum).variable(varNum).name():
                            cellMLRef['q' + str(speciesNum)] = []
                            cellMLRef['q' + str(speciesNum)].append(
                                newModels[modelName].component(compNum).variable(varNum).units())
                            cellMLRef['q' + str(speciesNum)].append(
                                newModels[modelName].component(compNum).variable(varNum).initialValue())
                            cellMLRef['q' + str(speciesNum)].append(
                                newModels[modelName].component(compNum).variable(varNum).id())

                if (modelName + '_' + newModels[modelName].component(compNum).variable(varNum).name() not in cellMLRef):
                    cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(varNum).name()] = []
                    cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(varNum).name()].append(
                        newModels[modelName].component(compNum).variable(varNum).units())
                    cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(varNum).name()].append(
                        newModels[modelName].component(compNum).variable(varNum).initialValue())
                    cellMLRef[modelName + '_' + newModels[modelName].component(compNum).variable(varNum).name()].append(
                        newModels[modelName].component(compNum).variable(varNum).id())

    unitNames = []
    for variable in cellMLRef:
        if list(cellMLRef[variable])[0].name() not in unitNames:
            unitNames.append(list(cellMLRef[variable])[0].name())
    common = []
    modelUnits = []
    for unit in unitNames:
        for variable in cellMLRef:
            if list(cellMLRef[variable])[0].name() == unit and list(cellMLRef[variable])[0].name() not in common:
                modelUnits.append(list(cellMLRef[variable])[0])
                common.append(list(cellMLRef[variable])[0].name())

    var_q = []
    for modelName in newModels:
        for x in [x for x in cellMLRef.keys() if modelName + '_' in x]:
            if modelName + '_' in x and cellMLRef[x][-1] in savedAnnotations[modelName].keys() and cellMLRef[x][
                0].name() == 'concentration_unit':
                # print(modelName,x,cellMLRef[x][-1])
                for species in speciesNoDuplicate:
                    if all(el in savedAnnotations[modelName][cellMLRef[x][-1]] for el in species[:-1]):
                        var_q.append((x, species[-1]))

    mathmlEquations = {}

    for modelName, model in newModels.items():
        for flux in V:
            if flux[0] == modelName and flux[1][0][1] == 'https://identifiers.org/opb/OPB_00592':  # is flux
                fluxID = flux[2]
        component = model.component(0)
        X1 = component.math().replace(
            '<math xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/2.0#">', '')
        X2 = X1.replace('<math xmlns="http://www.w3.org/1998/Math/MathML">', '')
        X3 = X2.replace('</math>', '')
        X4 = X3.replace('\n', '')
        X5 = X4.replace(' ', '')
        X6 = X5.replace('<apply><eq/><ci>' + fluxID + '</ci>',
                        '<apply><eq/><ci>' + modelName + '_' + fluxID + '</ci>')  # FIND IT BY CHECKING THE ANNOTATIONS
        X = X6.replace('cncellml', 'cn cellml')

        for compNum in range(newModels[modelName].componentCount()):
            for varNum in range(newModels[modelName].component(compNum).variableCount()):
                X = X.replace('>' + newModels[modelName].component(compNum).variable(varNum).name() + '<',
                              '>' + modelName + '_' + newModels[modelName].component(compNum).variable(
                                  varNum).name() + '<')

        mathmlEquations[modelName] = X

    preMath = ''
    newMath = ''
    for file_name in newModels:
        preMath += mathmlEquations[file_name]

    for speciesNum in odeRightSide:  # equivalent to speciesNoDuplicate
        newMath += q_to_mathml(odeRightSide[speciesNum], speciesNum, ode=True)
    mathBody = preMath + newMath

    model = libcellml.Model()
    model.setName("CompositeModel")

    for unit in modelUnits:
        model.addUnits(unit)

    # Create a new component
    component = libcellml.Component()
    component.setName("MyComponent")

    # Add the component to the model
    model.addComponent(component)

    v = libcellml.Variable()
    v.setName('t')
    v.setUnits('second')
    v.setInterfaceType("public")
    component.addVariable(v)

    # Create variables
    used = []
    for variable in cellMLRef:
        if variable not in [modelName + '_t' for modelName in newModels.keys()] and variable not in [modelName + '_X'
                                                                                                     for modelName in
                                                                                                     newModels.keys()]:
            if variable not in [x[0] for x in var_q]:
                v = libcellml.Variable()
                v.setName(variable)
                v.setUnits(list(cellMLRef[variable])[0])
                v.setInitialValue(list(cellMLRef[variable])[1])
                v.setInterfaceType("public")
                component.addVariable(v)
            else:
                for var, q in var_q:
                    if var == variable and q not in used:
                        v = libcellml.Variable()
                        v.setName(q)
                        v.setUnits(list(cellMLRef[variable])[0])
                        v.setInitialValue(list(cellMLRef[variable])[1])
                        v.setInterfaceType("public")
                        component.addVariable(v)
                        used.append(q)

    # Create the equation
    math_header = '<math xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/2.0#">'
    math_footer = '</math>'

    for modelName, modelinfo in newModels.items():
        mathBody = mathBody.replace('<ci>' + modelName + '_' + 't' + '</ci>', '<ci>' + 't' + '</ci>')
        mathBody = mathBody.replace(
            '<apply><eq/><apply><diff/><bvar><ci>t</ci></bvar><ci>' + modelName + '_' + 'X</ci></apply><cn cellml:units="dimensionless">0</cn></apply>',
            '')

    for var, q in var_q:
        mathBody = mathBody.replace('<ci>' + var + '</ci>', '<ci>' + q + '</ci>')

    component.setMath(math_header)
    component.appendMath(mathBody)
    component.appendMath(math_footer)

    annotator = Annotator()
    annotator.setModel(model)

    annotator.clearAllIds()
    annotator.assignAllIds()

    # Print the model
    printer = libcellml.Printer()
    # print(printer.printModel(model))

    writeCellML(model, printer, 'CompositeModel')

    return





def graphPlot(newCellmlNames, Species, speciesNoDuplicate, figSize):
    entities = [('H+', 'CHEBI:15378'), ('HCO3-', 'CHEBI:17544'), ('K+', 'CHEBI:29103'), ('Na+', 'CHEBI:29101'),
                ('Mg2+', 'CHEBI:18420'), ('Cl-', 'CHEBI:17996'), ('Ca2+', 'CHEBI:29108'), ('Fe2+', 'CHEBI:29033'),
                ('P', 'CHEBI:30207')]
    locations = [('Extracellular environment', 'fma70022'), ('Cytosol of stem cell', 'fma260697'),
                 ('Cytosol of neuron', 'fma226054'), ('Cytosol of connective tissue cell', 'fma260699'),
                 ('Cytosol of skeletal muscle fiber', 'fma226052'), ('Cytosol of cardiac myocyte', 'fma226050'),
                 ('Cytosol of striated muscle cell', 'fma226048'), ('Cytosol of smooth muscle cell', 'fma226046'),
                 ('Cytosol of muscle cell', 'fma226044'), ('Cytosol of hemal cell', 'fma260695'),
                 ('Cytosol of epithelial cell', 'fma260691')]

    connections = []
    SourcesTop = []
    SinksBelow = []
    for modelName in newCellmlNames:
        for contentSpecies in Species:
            if contentSpecies[0] == modelName:
                for x in contentSpecies[1]:
                    if 'isPartOf' in x[0]:
                        if 'fma70022' in x[1]:
                            SourcesTop.append([el[0] for el in entities if el[1] in [y[1] for y in contentSpecies[1] if
                                                                                     y[
                                                                                         0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and (
                                                                                                 'https://identifiers.org/chebi' in
                                                                                                 y[
                                                                                                     1] or 'https://identifiers.org/GO' in
                                                                                                 y[1])][0]][
                                                  0] + ' in ' + [el[0] for el in locations if el[1] in
                                                                 [y[1] for y in contentSpecies[1] if y[
                                                                     0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                              'http://identifiers.org/fma' in y[
                                                                          1] or 'https://identifiers.org/GO' in y[1])][
                                                                     0]][0])
                        else:
                            SinksBelow.append([el[0] for el in entities if el[1] in [y[1] for y in contentSpecies[1] if
                                                                                     y[
                                                                                         0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and (
                                                                                                 'https://identifiers.org/chebi' in
                                                                                                 y[
                                                                                                     1] or 'https://identifiers.org/GO' in
                                                                                                 y[1])][0]][
                                                  0] + ' in ' + [el[0] for el in locations if el[1] in
                                                                 [y[1] for y in contentSpecies[1] if y[
                                                                     0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                              'http://identifiers.org/fma' in y[
                                                                          1] or 'https://identifiers.org/GO' in y[1])][
                                                                     0]][0])

                    if x[1] == 'http://bime.uv.edu/semsim/Source':
                        A_disp = [el[0] for el in entities if el[1] in [y[1] for y in contentSpecies[1] if y[
                            0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and (
                                                                                    'https://identifiers.org/chebi' in
                                                                                    y[
                                                                                        1] or 'https://identifiers.org/GO' in
                                                                                    y[1])][0]][0] + ' in ' + \
                                 [el[0] for el in locations if el[1] in [y[1] for y in contentSpecies[1] if y[
                                     0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                                     'http://identifiers.org/fma' in y[
                                                                                 1] or 'https://identifiers.org/GO' in
                                                                                     y[1])][0]][0]
                        connections.append((A_disp, modelName))

                    if x[1] == 'http://bime.uv.edu/semsim/Sink':
                        B_disp = [el[0] for el in entities if el[1] in [y[1] for y in contentSpecies[1] if y[
                            0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and (
                                                                                    'https://identifiers.org/chebi' in
                                                                                    y[
                                                                                        1] or 'https://identifiers.org/GO' in
                                                                                    y[1])][0]][0] + ' in ' + \
                                 [el[0] for el in locations if el[1] in [y[1] for y in contentSpecies[1] if y[
                                     0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                                     'http://identifiers.org/fma' in y[
                                                                                 1] or 'https://identifiers.org/GO' in
                                                                                     y[1])][0]][0]
                        connections.append((modelName, B_disp))

    speciesNoDuplicateText = []
    for species in speciesNoDuplicate:
        speciesNoDuplicateText.append([el[0] for el in entities if el[1] in [y[1] for y in species if y[
            0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and ('https://identifiers.org/chebi' in y[
            1] or 'https://identifiers.org/GO' in y[1])][0]][0] + ' in ' + [el[0] for el in locations if el[1] in
                                                                            [y[1] for y in species if y[
                                                                                0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                                         'http://identifiers.org/fma' in
                                                                                         y[
                                                                                             1] or 'https://identifiers.org/GO' in
                                                                                         y[1])][0]][0])

    # Create an empty graph
    G = nx.DiGraph()

    # Draw the graph

    plt.figure(figsize=(figSize, figSize))

    fixed_positions = {key: ((figSize / len(newCellmlNames)) + i * (figSize / len(newCellmlNames)), 0) for key, i in
                       zip(newCellmlNames, range(len(newCellmlNames)))}
    source_positions = {key: ((figSize / len(SourcesTop)) + i * (figSize / len(SourcesTop)), 3) for key, i in
                        zip(SourcesTop, range(len(SourcesTop)))}
    sink_positions = {key: ((figSize / len(SinksBelow)) + i * (figSize / len(SinksBelow)), -3) for key, i in
                      zip(SinksBelow, range(len(SinksBelow)))}

    # Combine the positions
    positions = {**fixed_positions, **source_positions, **sink_positions}

    compartments = {}
    for species in speciesNoDuplicate:
        compartment = [el[0] for el in locations if el[1] in [y[1] for y in species if y[
            0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and ('http://identifiers.org/fma' in y[
            1] or 'https://identifiers.org/GO' in y[1])][0]][0]
        if compartment not in compartments.keys():
            compartments[compartment] = []
            compartments[compartment].append([el[0] for el in entities if el[1] in [y[1] for y in species if y[
                0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and ('https://identifiers.org/chebi' in y[
                1] or 'https://identifiers.org/GO' in y[1])][0]][0] + ' in ' + [el[0] for el in locations if el[1] in
                                                                                [y[1] for y in species if y[
                                                                                    0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                                             'http://identifiers.org/fma' in
                                                                                             y[
                                                                                                 1] or 'https://identifiers.org/GO' in
                                                                                             y[1])][0]][0])
        else:
            compartments[compartment].append([el[0] for el in entities if el[1] in [y[1] for y in species if y[
                0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and ('https://identifiers.org/chebi' in y[
                1] or 'https://identifiers.org/GO' in y[1])][0]][0] + ' in ' + [el[0] for el in locations if el[1] in
                                                                                [y[1] for y in species if y[
                                                                                    0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and (
                                                                                             'http://identifiers.org/fma' in
                                                                                             y[
                                                                                                 1] or 'https://identifiers.org/GO' in
                                                                                             y[1])][0]][0])

    # Create a color map with a unique color for each compartment
    compartment_colors = {compartment: color for compartment, color in
                          zip(compartments.keys(), cm.tab20(range(len(compartments.keys()))))}

    # Add nodes
    for compartment, nodes in compartments.items():
        G.add_nodes_from(nodes, compartment=compartment)
    G.add_nodes_from(newCellmlNames, compartment=compartment)

    # Draw the graph with different colors for compartments
    # for compartment,nodes in compartments.items():
    #     nx.draw_networkx_nodes(G, positions, nodelist=nodes, node_shape='o', node_color=compartment_colors[compartment], node_size=3000)
    for node, node_attrs in G.nodes(data=True):
        compartment = node_attrs['compartment']
        nx.draw_networkx_nodes(G, positions, nodelist=[node],
                               node_shape='o', node_color=compartment_colors[compartment], node_size=3000)

    # Draw transporters (model names)
    nx.draw_networkx_nodes(G, positions, nodelist=newCellmlNames, node_shape='H', node_color='lightgrey',
                           node_size=5000)

    # Add edges
    G.add_edges_from(connections)

    edge_arrows = nx.draw_networkx_edges(G, positions, edge_color='grey', arrowsize=30, arrowstyle="-|>", width=0.5,
                                         arrows=True)  # ,  connectionstyle='angle,angleA=0,angleB=90')

    for arrow in edge_arrows:
        arrow.set_zorder(2)

    labels = nx.draw_networkx_labels(G, positions, font_color='black', font_size=10)  # , font_weight='bold')

    # Adjust the positions of labels to avoid overlap
    adjust_text(texts=list(labels.values()), autoalign='xy', va='center')

    # Show the graph
    plt.show()





def templateGraphPlot(Type, figSize, equations):
    types = ['Template 3', 'Template 4', 'Template 7', 'Template 11', 'Template 5']

    elements = dict((el, []) for el in types)
    elements['Template 3'] = ['$q_{A_o}$', '$q_{A_i}$']
    elements['Template 4'] = ['$q_{A_o}$', '$q_{A_i}$', '$q_{B_o}$', '$q_{B_i}$']
    elements['Template 7'] = ['$q_{A_o}$', '$q_{A_i}$', '$q_{B_o}$', '$q_{B_i}$', '$q_{C_o}$', '$q_{C_i}$']
    elements['Template 11'] = ['$q_{A_o}$', '$q_{A_i}$', '$q_{B_o}$', '$q_{B_i}$', '$q_{C_o}$', '$q_{C_i}$',
                               '$q_{D_o}$', '$q_{D_i}$']
    elements['Template 5'] = ['$q_{A_o}$', '$q_{A_i}$', '$q_{B_o}$', '$q_{B_i}$', '$V_m$']

    connections = dict((el, []) for el in types)
    connections['Template 3'] = [('$q_{A_o}$', 'Template 3'), ('Template 3', '$q_{A_i}$')]
    connections['Template 4'] = [('$q_{A_o}$', 'Template 4'), ('Template 4', '$q_{A_i}$'), ('$q_{B_i}$', 'Template 4'),
                                 ('Template 4', '$q_{B_o}$')]
    connections['Template 7'] = [('$q_{A_i}$', 'Template 7'), ('Template 7', '$q_{A_o}$'), ('$q_{B_o}$', 'Template 7'),
                                 ('Template 7', '$q_{B_i}$'), ('$q_{C_o}$', 'Template 7'), ('Template 7', '$q_{C_i}$')]
    connections['Template 11'] = [('$q_{A_i}$', 'Template 11'), ('Template 11', '$q_{A_o}$'),
                                  ('$q_{B_o}$', 'Template 11'), ('Template 11', '$q_{B_i}$'),
                                  ('$q_{C_o}$', 'Template 11'), ('Template 11', '$q_{C_i}$'),
                                  ('$q_{D_o}$', 'Template 11'), ('Template 11', '$q_{D_i}$')]
    connections['Template 5'] = [('$q_{A_o}$', 'Template 5'), ('Template 5', '$q_{A_i}$'), ('$q_{B_o}$', 'Template 5'),
                                 ('Template 5', '$q_{B_i}$'), ('$V_m$', 'Template 5')]

    SourcesTop = []
    SinksBelow = []
    membrane = []

    for x in elements[Type]:
        if x in ['$q_{A_i}$', '$q_{B_i}$', '$q_{C_i}$', '$q_{D_i}$']:
            SinksBelow.append(x)
            # connections.append((Type,x))
        elif x == '$V_m$':
            membrane.append(x)
            # connections.append((x,Type))
        else:
            SourcesTop.append(x)
            # connections.append((x,Type))

    # Create an empty graph
    G = nx.DiGraph()

    # Add nodes

    G.add_nodes_from(elements[Type])
    G.add_nodes_from([Type])
    try:
        if membrane:
            G.add_nodes_from(membrane)
    except:
        pass

    # Add edges
    G.add_edges_from(connections[Type])

    # Draw the graph

    plt.figure(figsize=(figSize, figSize))

    fixed_positions = {key: ((figSize / len([Type])) + i * (figSize / len([Type])), 0) for key, i in
                       zip([Type], range(len([Type])))}
    source_positions = {key: ((figSize / len(SourcesTop + membrane)) + i * (figSize / len(SourcesTop + membrane)), 3)
                        for key, i in zip(SourcesTop + membrane, range(len(SourcesTop + membrane)))}
    sink_positions = {key: ((figSize / len(SinksBelow)) + i * (figSize / len(SinksBelow)), -3) for key, i in
                      zip(SinksBelow, range(len(SinksBelow)))}

    # Combine the positions
    positions = {**fixed_positions, **source_positions, **sink_positions}

    # Draw nodes with custom shapes
    for node in [Type]:
        nx.draw_networkx_nodes(G, positions, nodelist=[node], node_shape='H', node_color='orange', node_size=5000)
    for node in SourcesTop + SinksBelow:
        nx.draw_networkx_nodes(G, positions, nodelist=[node], node_shape='o', node_color='lightblue', node_size=3000)
    for node in membrane:
        nx.draw_networkx_nodes(G, positions, nodelist=[node], node_shape='v', node_color='yellowgreen', node_size=3000)

    edge_arrows = nx.draw_networkx_edges(G, positions, edge_color='gray', arrowsize=30, arrowstyle="-|>", width=0.5,
                                         arrows=True)

    for arrow in edge_arrows:
        arrow.set_zorder(2)

    labels = nx.draw_networkx_labels(G, positions, font_color='black', font_size=10)

    # Adjust the positions of labels to avoid overlap
    adjust_text(texts=list(labels.values()), autoalign='y', va='center')

    # Show the graph
    plt.show()





def annotate(key, element, selectedAnnotations):
    # Create a dropdown widget for selecting the search column
    column_names_left = [('H+', 'entity$CHEBI:15378'), ('HCO3-', 'entity$CHEBI:17544'), ('K+', 'entity$CHEBI:29103'),
                         ('Na+', 'entity$CHEBI:29101'), ('Mg2+', 'entity$CHEBI:18420'), ('Cl-', 'entity$CHEBI:17996'),
                         ('Ca2+', 'entity$CHEBI:29108'), ('Fe2+', 'entity$CHEBI:29033'), ('P', 'entity$CHEBI:30207')]
    column_names_widget_left = widgets.Dropdown(
        options=column_names_left,
        style={'description_width': 'initial'},
        disabled=False,
        layout=Layout(width='400px', height='60px')
    )

    column_names_middle = [('Extracellular environment', 'isPartOf$fma70022'),
                           ('Cytosol of stem cell', 'isPartOf$fma260697'), ('Cytosol of neuron', 'isPartOf$fma226054'),
                           ('Cytosol of connective tissue cell', 'isPartOf$fma260699'),
                           ('Cytosol of skeletal muscle fiber', 'isPartOf$fma226052'),
                           ('Cytosol of cardiac myocyte', 'isPartOf$fma226050'),
                           ('Cytosol of striated muscle cell', 'isPartOf$fma226048'),
                           ('Cytosol of smooth muscle cell', 'isPartOf$fma226046'),
                           ('Cytosol of muscle cell', 'isPartOf$fma226044'),
                           ('Cytosol of hemal cell', 'isPartOf$fma260695'),
                           ('Cytosol of epithelial cell', 'isPartOf$fma260691')]
    column_names_widget_middle = widgets.Dropdown(
        options=column_names_middle,
        style={'description_width': 'initial'},
        disabled=False,
        layout=Layout(width='400px', height='60px')
    )

    ## Set the font size using CSS
    font_size = '28px'  # Font size in CSS format

    # Generate the CSS code to update the font size
    css_code = f"""
    <style>
    .widget-dropdown select {{
        font-size: {font_size};
    }}
    .widget-dropdown .widget-label {{
        font-size: {font_size};
    }}
    </style>
    """

    def left_dropdown_eventhandler(change1):
        for i in selectedAnnotations[key][element]:
            if 'entity$' in i:
                selectedAnnotations[key][element].remove(i)
        selectedAnnotations[key][element].append(change1.new)
        with open('./temporary files/selectedAnnotations.json', 'w') as f1:
            json.dump(selectedAnnotations, f1, indent=4, sort_keys=True, ensure_ascii=False)

    def middle_dropdown_eventhandler(change2):
        for i in selectedAnnotations[key][element]:
            if 'isPartOf$' in i:
                selectedAnnotations[key][element].remove(i)
        selectedAnnotations[key][element].append(change2.new)
        with open('./temporary files/selectedAnnotations.json', 'w') as f1:
            json.dump(selectedAnnotations, f1, indent=4, sort_keys=True, ensure_ascii=False)

    column_names_widget_left.observe(left_dropdown_eventhandler, names='value')
    column_names_widget_middle.observe(middle_dropdown_eventhandler, names='value')

    # Display the widgets
    text_left = widgets.HTML(value="<h1><b><font color='salmon'>Physical entity:<b><h1>")
    display(text_left)
    display(column_names_widget_left)
    text_middle = widgets.HTML(value="<h1><b><font color='salmon'>Contained in:<b><h1>")
    display(text_middle)
    display(column_names_widget_middle)
    display(HTML(css_code))

    return selectedAnnotations





def calculate_v(pVals, *species):
    f = open("./temporary files/pVals.json")
    pVals = json.load(f)

    if singleSelection.value == 'Template 3':
        Ao, Ai = species
        return (pVals[singleSelection.value]['p1'] * Ao - pVals[singleSelection.value]['p2'] * Ai) / (
                    pVals[singleSelection.value]['p3'] + pVals[singleSelection.value]['p4'] * Ao +
                    pVals[singleSelection.value]['p5'] * Ai * Ao + pVals[singleSelection.value]['p6'] * Ai)
    if singleSelection.value == 'Template 4':
        Ao, Ai, Bo, Bi = species
        return (pVals[singleSelection.value]['p1'] * Ao * Bi - pVals[singleSelection.value]['p2'] * Ai * Bo) / (
                    pVals[singleSelection.value]['p3'] * Ao + pVals[singleSelection.value]['p4'] * Ai +
                    pVals[singleSelection.value]['p5'] * Bo + pVals[singleSelection.value]['p6'] * Bi +
                    pVals[singleSelection.value]['p7'] * Ao * Bi + pVals[singleSelection.value]['p8'] * Ai * Ao +
                    pVals[singleSelection.value]['p9'] * Ai * Bo + pVals[singleSelection.value]['p10'] * Bo * Bi)
    if singleSelection.value == 'Template 7':
        Ao, Ai, Bo, Bi, Co, Ci = species
        return (pVals[singleSelection.value]['p1'] * Ai * Bo * Co - pVals[singleSelection.value][
            'p2'] * Ao * Bi * Ci) / (pVals[singleSelection.value]['p3'] * Ai + pVals[singleSelection.value]['p4'] * Ao +
                                     pVals[singleSelection.value]['p5'] * Bi + pVals[singleSelection.value]['p6'] * Bo +
                                     pVals[singleSelection.value]['p7'] * Ci + pVals[singleSelection.value]['p8'] * Co +
                                     pVals[singleSelection.value]['p9'] * Ai * Ao + pVals[singleSelection.value][
                                         'p10'] * Bi * Bo + pVals[singleSelection.value]['p11'] * Ci * Co +
                                     pVals[singleSelection.value]['p12'] * Ai * Bo * Co + pVals[singleSelection.value][
                                         'p13'] * Bi * Ci * Ao)
    if singleSelection.value == 'Template 11':
        Ao, Ai, Bo, Bi, Co, Ci, Do, Di = species
        return (pVals[singleSelection.value]['p1'] * Ai * Bo * Co * Do - pVals[singleSelection.value][
            'p2'] * Ao * Bi * Ci * Di) / (
                           pVals[singleSelection.value]['p3'] * Ai + pVals[singleSelection.value]['p4'] * Ao +
                           pVals[singleSelection.value]['p5'] * Bi + pVals[singleSelection.value]['p6'] * Bo +
                           pVals[singleSelection.value]['p7'] * Ci + pVals[singleSelection.value]['p8'] * Co +
                           pVals[singleSelection.value]['p9'] * Di + pVals[singleSelection.value]['p10'] * Do +
                           pVals[singleSelection.value]['p11'] * Ai * Ao + pVals[singleSelection.value][
                               'p12'] * Bi * Bo + pVals[singleSelection.value]['p13'] * Ci * Co +
                           pVals[singleSelection.value]['p14'] * Di * Do + pVals[singleSelection.value][
                               'p15'] * Ai * Bo * Co * Do + pVals[singleSelection.value]['p16'] * Bi * Ci * Ao * Di)
    if singleSelection.value == 'Template 5':
        Ao, Ai, Bo, Bi, V_m = species
        F = 96485
        R = 8.31
        T = 310
        return (pVals[singleSelection.value]['p1'] * Bo * Ao * exp(F * V_m / (R * T)) - pVals[singleSelection.value][
            'p2'] * Ai * Bi) / (pVals[singleSelection.value]['p11'] * Bi * Ai * Ao * exp(F * V_m / (R * T)) +
                                pVals[singleSelection.value]['p12'] * Bi * Bo * Ai * Ao * exp(F * V_m / (R * T)) +
                                pVals[singleSelection.value]['p3'] * Bi * Bo + pVals[singleSelection.value][
                                    'p4'] * Bo * Ai + pVals[singleSelection.value]['p5'] * Bo * Ao * exp(
                    F * V_m / (R * T)) + pVals[singleSelection.value]['p6'] * Bi * Ai + pVals[singleSelection.value][
                                    'p7'] * Bi * Ao * exp(F * V_m / (R * T)) + pVals[singleSelection.value][
                                    'p8'] * Bi * Bo * Ai + pVals[singleSelection.value]['p9'] * Bi * Bo * Ao * exp(
                    F * V_m / (R * T)) + pVals[singleSelection.value]['p10'] * Bo * Ai * Ao * exp(F * V_m / (R * T)))





def update_figure(control_variable, control_value, stepTime, timespan, **const_values):
    time = np.linspace(0, timespan, 1000)
    f = open("./temporary files/pVals.json")
    pVals = json.load(f)

    V = []
    if singleSelection.value == 'Template 3':
        if control_variable == 'Ao':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ao']
                    else:
                        return control_value

                V.append(calculate_v(pVals, step_func(t), const_values['Ai']))

        if control_variable == 'Ai':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ai']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], step_func(t)))

    if singleSelection.value == 'Template 4':
        if control_variable == 'Ao':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ao']
                    else:
                        return control_value

                V.append(calculate_v(pVals, step_func(t), const_values['Ai'], const_values['Bo'], const_values['Bi']))

        if control_variable == 'Ai':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ai']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], step_func(t), const_values['Bo'], const_values['Bi']))

        if control_variable == 'Bo':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bo']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], step_func(t), const_values['Bi']))

        if control_variable == 'Bi':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bi']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], step_func(t)))

    if singleSelection.value == 'Template 7':
        if control_variable == 'Ao':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ao']
                    else:
                        return control_value

                V.append(calculate_v(pVals, step_func(t), const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                     const_values['Co'], const_values['Ci']))

        if control_variable == 'Ai':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ai']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], step_func(t), const_values['Bo'], const_values['Bi'],
                                     const_values['Co'], const_values['Ci']))

        if control_variable == 'Bo':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bo']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], step_func(t), const_values['Bi'],
                                     const_values['Co'], const_values['Ci']))

        if control_variable == 'Bi':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bi']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], step_func(t),
                                     const_values['Co'], const_values['Ci']))

        if control_variable == 'Co':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Co']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                step_func(t), const_values['Ci']))

        if control_variable == 'Ci':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ci']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                const_values['Co'], step_func(t)))

    if singleSelection.value == 'Template 11':
        if control_variable == 'Ao':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ao']
                    else:
                        return control_value

                V.append(calculate_v(pVals, step_func(t), const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                     const_values['Co'], const_values['Ci'], const_values['Do'], const_values['Di']))

        if control_variable == 'Ai':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ai']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], step_func(t), const_values['Bo'], const_values['Bi'],
                                     const_values['Co'], const_values['Ci'], const_values['Do'], const_values['Di']))

        if control_variable == 'Bo':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bo']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], step_func(t), const_values['Bi'],
                                     const_values['Co'], const_values['Ci'], const_values['Do'], const_values['Di']))

        if control_variable == 'Bi':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bi']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], step_func(t),
                                     const_values['Co'], const_values['Ci'], const_values['Do'], const_values['Di']))

        if control_variable == 'Co':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Co']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                step_func(t), const_values['Ci'], const_values['Do'], const_values['Di']))

        if control_variable == 'Ci':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ci']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                const_values['Co'], step_func(t), const_values['Do'], const_values['Di']))

        if control_variable == 'Do':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Do']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                const_values['Co'], const_values['Ci'], step_func(t), const_values['Di']))

        if control_variable == 'Di':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Di']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                const_values['Co'], const_values['Ci'], const_values['Do'], step_func(t)))

    if singleSelection.value == 'Template 5':
        if control_variable == 'Ao':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ao']
                    else:
                        return control_value

                V.append(calculate_v(pVals, step_func(t), const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                     const_values['V_m']))

        if control_variable == 'Ai':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Ai']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], step_func(t), const_values['Bo'], const_values['Bi'],
                                     const_values['V_m']))

        if control_variable == 'Bo':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bo']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], step_func(t), const_values['Bi'],
                                     const_values['V_m']))

        if control_variable == 'Bi':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['Bi']
                    else:
                        return control_value

                V.append(calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], step_func(t),
                                     const_values['V_m']))

        if control_variable == 'V_m':
            for t in time:
                def step_func(t):
                    if t < stepTime:
                        return const_values['V_m']
                    else:
                        return control_value

                V.append(
                    calculate_v(pVals, const_values['Ao'], const_values['Ai'], const_values['Bo'], const_values['Bi'],
                                step_func(t)))

    # Update the plot
    matplotlib.rcParams.update({'font.size': 15})
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(time, V, color='m', linewidth=5, label='$V_{ss}$')
    ax.set_xlabel('Time (s)')
    ax.grid(color='gray', linestyle='-', linewidth=0.5)
    ax.set_ylabel('Flux ($u$mol/s)')
    ax.legend(loc='best')
    formatted_control_value = f'{control_value:.2f}'
    formatted_control_variable = [f'$q_{control_variable}$' if control_variable != 'V_m' else f'${control_variable}$'][
        0]
    ax.set_title(f'Control variable: {formatted_control_variable}, Control value: {formatted_control_value}')
    plt.show()




#
# <div class="alert alert-block alert-warning">
#     <h1><b>Model Construction</b></h1>
# </div>
#
# ***




savedModels = []
savedAnnotations = {}
savedkfkr = {}
savedtypes = {}
savedSpeciesNoDuplicateSingleDict = {}
initialValues = {}
Nmatrices = {}
savedpVals = {}

equations = {}
equations[
    'Template 3'] = r"$V_{SS} = \frac{(p_{1} q_{Ao} - p_{2} q_{Ai})}{(p_{3} + p_{4} q_{Ao} + p_{5} q_{Ai} q_{Ao} + p_{6} q_{Ai})}$"
equations[
    'Template 4'] = r"$V_{SS} = \frac{(p_{1} q_{Ao} q_{Bi} - p_{2} q_{Ai} q_{Bo})}{(p_{3} q_{Ao} + p_{4} q_{Ai} + p_{5} q_{Bo} + p_{6} q_{Bi} + p_{7} q_{Ao} q_{Bi} + p_{8} q_{Ai} q_{Ao} + p_{9} q_{Ai} q_{Bo} + p_{10} q_{Bo} q_{Bi})}$"
equations[
    'Template 7'] = r"$V_{SS} = \frac{(p_{1} q_{Ai} q_{Bo} q_{Co} - p_{2} q_{Ao} q_{Bi} q_{Ci})}{(p_{3} q_{Ai} + p_{4} q_{Ao} + p_{5} q_{Bi} + p_{6} q_{Bo} + p_{7} q_{Ci} + p_{8} q_{Co} + p_{9} q_{Ai} q_{Ao} + p_{10} q_{Bi} q_{Bo} + p_{11} q_{Ci} q_{Co} + p_{12} q_{Ai} q_{Bo} q_{Co} + p_{13} q_{Bi} q_{Ci} q_{Ao} )}$"
equations[
    'Template 11'] = r"$V_{SS} = \frac{(p_{1} q_{Ai} q_{Bo} q_{Co} q_{Do} - p_{2} q_{Ao} q_{Bi} q_{Ci} q_{Di})}{(p_{3} q_{Ai} + p_{4} q_{Ao} + p_{5} q_{Bi} + p_{6} q_{Bo} + p_{7} q_{Ci} + p_{8} q_{Co} + p_{9} q_{Di} + p_{10} q_{Do} + p{11} q_{Ai} q_{Ao} + p_{12} q_{Bi} q_{Bo} + p_{13} q_{Ci} q_{Co} + p_{14} q_{Di} q_{Do} + p_{15} q_{Ai} q_{Bo} q_{Co} q_{Do} + p_{16} q_{Bi} q_{Ci} q_{Ao} q_{Di})}$"
equations[
    'Template 5'] = r"$V_{SS} = \frac{(p_{1} q_{Bo} q_{Ao} e^{(F V_m/R T)}-p_{2} q_{Ai} q_{Bi})}{(p_{3} q_{Bi} q_{Bo}+p_{4} q_{Bo} q_{Ai} +p_{5} q_{Bo} q_{Ao} e^{(F V_m/R T)}+p_{6} q_{Bi} q_{Ai}+p_{7} q_{Bi} q_{Ao} e^{(F V_m/R T)}+p_{8} q_{Bi} q_{Bo} q_{Ai}+p_{9} q_{Bi} q_{Bo} q_{Ao} e^{(F V_m/R T)}+p_{10} q_{Bo} q_{Ai} q_{Ao} e^{(F V_m/R T)} +p_{11} q_{Bi} q_{Ai} q_{Ao} e^{(F V_m/R T)}+p_{12} q_{Bi} q_{Bo} q_{Ai} q_{Ao} e^{(F V_m/R T)})}$"

# # Section 1: Template selection




# Template Selection and BG structure demonstration

singleSelection = widgets.Dropdown(
    options=['Template 3', 'Template 4', 'Template 7', 'Template 11', 'Template 5'],
    description='Select transporter type:',
    style={'description_width': 'initial', 'description_height': 'initial'},
    disabled=False,
    layout=Layout(width='700px', height='60px')
)

## Set the font size using CSS
font_size = '28px'  # Font size in CSS format

# Generate the CSS code to update the font size
custom_css = f"""
<style>
.custom-text-style .widget-label {{
    font-size: {font_size};
    height: initial;
    width: initial;
}}    
.widget-dropdown select {{
    font-size: {font_size};
}}
.widget-dropdown .widget-label {{
    font-size: {font_size};
}}

</style>
"""

button_load = widgets.Button(
    description='Load template',
    tooltip='Load',
    style={'description_width': 'initial'},
    button_style='info',
    layout=Layout(width='300px', height='60px')
)
button_load.style.font_family = 'monospace'
button_load.style.font_size = '28px'
output1 = widgets.Output(layout={'border': '1px solid gray'})
outputMathTemp = widgets.Output(layout={'border': '1px solid gray'})


def button_load_clicked(event):
    with output1:
        clear_output()
        templateGraphPlot(singleSelection.value, 10, equations)
    with outputMathTemp:
        clear_output()

        # Create a plot to display the expression
        fig, ax = plt.subplots(figsize=(6, 1))
        ax.text(0.5, 0.5, equations[singleSelection.value], fontsize=20, ha='center', va='center')
        ax.axis('off')
        plt.show()

        # # Create a Math object with the expression
        # math_expression = Math(equations[singleSelection.value])
        # # Display the Math object
        # display(math_expression)


# Apply the custom CSS styling
singleSelection.add_class('widget-dropdown')
singleSelection.add_class('custom-text-style')

button_load.on_click(button_load_clicked)
vbox_result = widgets.VBox([button_load, output1, outputMathTemp])

vbox_text = widgets.VBox([singleSelection, vbox_result])
page1 = widgets.HBox([vbox_text])

display(page1)
display(HTML(custom_css))

# # Section 2: Model Description




types = ['Template 3', 'Template 4', 'Template 7', 'Template 11', 'Template 5']

valueRequired1 = dict((el, []) for el in types)
valueRequired1['Template 3'] = ['Ao', 'Ai']
valueRequired1['Template 4'] = ['Ao', 'Ai', 'Bo', 'Bi']
valueRequired1['Template 7'] = ['Ao', 'Ai', 'Bo', 'Bi', 'Co', 'Ci']
valueRequired1['Template 11'] = ['Ao', 'Ai', 'Bo', 'Bi', 'Co', 'Ci', 'Do', 'Di']
valueRequired1['Template 5'] = ['Ao', 'Ai', 'Bo', 'Bi']

button_annotate = widgets.Button(
    description='Describe the model!',
    tooltip='Description',
    style={'description_width': 'initial'},
    button_style='info',
    layout=Layout(width='400px', height='60px')
)
button_annotate.style.font_family = 'monospace'
button_annotate.style.font_size = '28px'
output_annotations = widgets.Output(layout={'height': '100%'})


def button_annotate_clicked(event):
    with output_annotations:
        output_annotations.clear_output()
        text00 = widgets.HTML(
            value="<h3><b><font color='maroon'>Annotate the following components as in the loaded template:<b><h3>")
        display(text00)
        selectedAnnotations = {}

        key = singleSelection.value
        selectedAnnotations[key] = {}

        for element in valueRequired1[singleSelection.value]:
            selectedAnnotations[key]['V_SS'] = ['isVersionOf$OPB_00592']
            selectedAnnotations[key][element] = ['isVersionOf$OPB_00425', 'entity$CHEBI:15378', 'isPartOf$fma70022']

            if key == 'Template 3':
                if element == 'Ao':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Ai':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
            if key == 'Template 4':
                if element == 'Ao':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Ai':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
                if element == 'Bi':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Bo':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
            if key == 'Template 7':
                if element == 'Ao':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
                if element == 'Ai':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Bo':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Bi':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
                if element == 'Co':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Ci':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
            if key == 'Template 11':
                if element == 'Ao':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
                if element == 'Ai':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Bo':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Bi':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
                if element == 'Co':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Ci':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
                if element == 'Do':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Di':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')
            if key == 'Template 5':
                if element == 'Ao':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Source')
                if element == 'Ai':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Sink')
                if element == 'Bo':
                    selectedAnnotations[key][element].append('hasSourceParticipant$Source')
                if element == 'Bi':
                    selectedAnnotations[key][element].append('hasSinkParticipant$Sink')

            text = widgets.HTML(value="<h1><b><font color='teal'>q<sub>{}</sub><b><h1>".format(element))

            display(text)

            selectedAnnotations = annotate(key, element, selectedAnnotations)

            # Save the values in .json file
        with open("./temporary files/selectedAnnotations.json", "w") as outfile:
            json.dump(selectedAnnotations, outfile, indent=4, sort_keys=True, ensure_ascii=False)


button_annotate.on_click(button_annotate_clicked)
annoat_page = widgets.VBox([button_annotate, output_annotations])
display(annoat_page)

# # Section 3: Parameterization




CellVol_button = widgets.Button(description='Display', layout=Layout(width='300px', height='60px'))
CellVol_button.style.font_size = '28px'
CellVol_button.style.button_color = 'skyblue'
CellVol_button.style.font_family = 'monospace'
CellVol_output = widgets.Output()

cellVol = widgets.FloatText(value=1, description='Insert the cell volume: ', style={'description_width': 'initial'},
                            layout=Layout(width='450px', height='60px'))
cellVol.style.font_family = 'monospace'
# cellUnit = widgets.Dropdown(options=[('femtolitre',1e-15), ('picolitre',1e-12), ('nanolitre',1e-9), ('microlitre',1e-6)],  style={'description_width': 'initial'})
cellUnit = widgets.Dropdown(options=['femtolitre', 'picolitre', 'nanolitre', 'microlitre'],
                            style={'description_width': 'initial'}, layout=Layout(width='250px', height='60px'))

vbox_result = widgets.HBox([cellVol, cellUnit])
display(vbox_result)


def validate_positive_value(change):
    value = change['new']
    if value < 0:
        change['owner'].value = 0


cellVol.observe(validate_positive_value, 'value')


def CellVol_button_clicked(event):
    with CellVol_output:
        clear_output()

        # formatted_cellVol = float("{:.2f}".format(cellVol.value))* cellUnit.value
        # print('The cell volume is: '+str(formatted_cellVol)+' (Litre)')
        print('The cell volume is: ' + str(cellVol.value) + ' ' + cellUnit.value)


# Apply the custom CSS class to the FloatText widget
cellVol.add_class('custom-text-style')

CellVol_button.on_click(CellVol_button_clicked)
vbox_result = widgets.VBox([CellVol_button, CellVol_output])
display(vbox_result)
display(HTML(custom_css))




# Getting the SS values
types = ['Template 3', 'Template 4', 'Template 7', 'Template 11', 'Template 5']
valueRequired = dict((el, []) for el in types)

valueRequired['Template 3'] = ['V_SS', 'Ao', 'Ai']
valueRequired['Template 4'] = ['V_SS', 'Ao', 'Ai', 'Bo', 'Bi']
valueRequired['Template 7'] = ['V_SS', 'Ao', 'Ai', 'Bo', 'Bi', 'Co', 'Ci']
valueRequired['Template 11'] = ['V_SS', 'Ao', 'Ai', 'Bo', 'Bi', 'Co', 'Ci', 'Do', 'Di']
valueRequired['Template 5'] = ['V_SS', 'Ao', 'Ai', 'Bo', 'Bi', 'V_m']

button_SS = widgets.Button(
    description='Steady-state Values',
    tooltip='Description',
    style={'description_width': 'initial'},
    button_style='info',
    layout=Layout(width='400px', height='60px')
)
button_SS.style.font_family = 'monospace'
button_SS.style.font_size = '28px'

buttonExperimentalData = widgets.Button(description='Experimental Data',
                                        tooltip='Description',
                                        style={'description_width': 'initial'},
                                        button_style='info',
                                        layout=Layout(width='400px', height='60px'))
buttonExperimentalData.style.font_family = 'monospace'
buttonExperimentalData.style.font_size = '28px'

output_parameterization = widgets.Output(layout={})
output_csv = widgets.Output(layout={})


# Tab 1 functionality
def button_SS_clicked(event):
    with output_parameterization:
        clear_output()

        file = open('./temporary files/selectedAnnotations.json')
        annots = json.load(file)

        entities = [('H+', 'CHEBI:15378'), ('HCO3-', 'CHEBI:17544'), ('K+', 'CHEBI:29103'), ('Na+', 'CHEBI:29101'),
                    ('Mg2+', 'CHEBI:18420'), ('Cl-', 'CHEBI:17996'), ('Ca2+', 'CHEBI:29108'), ('Fe2+', 'CHEBI:29033'),
                    ('P', 'CHEBI:30207')]
        locations = [('Extracellular environment', 'fma70022'), ('Cytosol of stem cell', 'fma260697'),
                     ('Cytosol of neuron', 'fma226054'), ('Cytosol of connective tissue cell', 'fma260699'),
                     ('Cytosol of skeletal muscle fiber', 'fma226052'), ('Cytosol of cardiac myocyte', 'fma226050'),
                     ('Cytosol of striated muscle cell', 'fma226048'), ('Cytosol of smooth muscle cell', 'fma226046'),
                     ('Cytosol of muscle cell', 'fma226044'), ('Cytosol of hemal cell', 'fma260695'),
                     ('Cytosol of epithelial cell', 'fma260691')]

        global floatSSget, Description
        floatSSget = [];
        Description = []
        for i, el in enumerate(valueRequired[singleSelection.value]):
            if el == 'V_SS':
                Description.append(widgets.HTML(value='Steady-state flux'))
            elif el == 'V_m':
                Description.append(widgets.HTML(value='Membrane voltage'))
            else:
                Description.append(widgets.HTML(value=[y[0] for y in entities if y[1] in
                                                       [x for x in annots[singleSelection.value][el] if 'entity$' in x][
                                                           0]][0] + ' in ' + [y[0] for y in locations if y[1] in
                                                                              [x for x in
                                                                               annots[singleSelection.value][el] if
                                                                               'isPartOf$' in x][0]][0]))
            floatSSget.append(widgets.FloatText(
                value=0.0,
                description=['$q_{}$'.format(el) if el not in ['V_SS', 'V_m'] else el][0] +
                            ['  ($u$mol/s)' if el == 'V_SS' else '  (mV)' if el == 'V_m' else '  ($u$mol)'][0],
                disabled=False,
                style={'description_width': 'initial'},
                layout=Layout(width='300px', height='60px')
            ))
            for F in floatSSget:
                F.add_class('widget-dropdown')
                F.add_class('custom-text-style')

        def validate_positive_value(change):
            value = change['new']
            if value < 0:
                change['owner'].value = 0.0

        for i, el in enumerate(valueRequired[singleSelection.value]):
            if el != 'V_m':
                floatSSget[i].observe(validate_positive_value, 'value')

        button_saveSS = widgets.Button(
            description='Save',
            tooltip='Save',
            style={'description_width': 'initial'},
            layout=Layout(width='120px', height='40px'))
        button_saveSS.style.font_family = 'monospace'
        button_saveSS.style.font_size = '24px'
        button_saveSS.style.button_color = 'lightgreen'

        def button_saveSS_clicked(event):
            with output_parameterization:
                # clear_output()
                text = widgets.HTML(value="<h4><b>The inserted steady-state values have been saved.<b><h4>")
                display(text)

                for floatSS in floatSSget:
                    if floatSS.description != 'V_m  (mV)':
                        floatSS.value = copy.copy(floatSS.value / cellVol.value)

                pVals = pValsFunc(singleSelection, floatSSget, dataset=False)
                # Replace the parameter values in the equation
                equation = copy.deepcopy(equations[singleSelection.value])
                for param, value in pVals[singleSelection.value].items():
                    param_pattern = re.escape('p_{' + param[1:] + '}')
                    if abs(value) >= 0.01:
                        formatted_value = "{:.2f}".format(value)
                    else:
                        formatted_value = "{:.2e}".format(value)
                    equation = re.sub(param_pattern, formatted_value, equation)

                # Create a plot to display the expression
                fig, ax = plt.subplots(figsize=(6, 1))
                ax.text(0.5, 0.5, equation, fontsize=20, ha='center', va='center')
                ax.axis('off')
                plt.show()

        button_saveSS.on_click(button_saveSS_clicked)

        for i in range(len(Description)):
            display(Description[i], floatSSget[i])
            display(HTML(custom_css))

        display(button_saveSS)
        display(HTML(custom_css))


# Tab 2 functionality
def buttonExperimentalData_clicked(event):
    with output_csv:
        clear_output()
        # Create widgets for the second tab
        file_upload_widget = widgets.FileUpload(accept='.csv')
        submit_button_widget = widgets.Button(description='Submit')

        def submit_button_clicked(event):
            with output_csv:
                clear_output()
                # Get the uploaded file
                uploaded_files = file_upload_widget.value[0]
                if uploaded_files:
                    file_contents = uploaded_files['content']
                    df = pd.read_csv(io.BytesIO(file_contents))

                    ###############################################################################################################

                    file = open('./temporary files/selectedAnnotations.json')
                    annots = json.load(file)

                    entities = [('H+', 'CHEBI:15378'), ('HCO3-', 'CHEBI:17544'), ('K+', 'CHEBI:29103'),
                                ('Na+', 'CHEBI:29101'), ('Mg2+', 'CHEBI:18420'), ('Cl-', 'CHEBI:17996'),
                                ('Ca2+', 'CHEBI:29108'), ('Fe2+', 'CHEBI:29033'), ('P', 'CHEBI:30207')]
                    locations = [('Extracellular environment', 'fma70022'), ('Cytosol of stem cell', 'fma260697'),
                                 ('Cytosol of neuron', 'fma226054'), ('Cytosol of connective tissue cell', 'fma260699'),
                                 ('Cytosol of skeletal muscle fiber', 'fma226052'),
                                 ('Cytosol of cardiac myocyte', 'fma226050'),
                                 ('Cytosol of striated muscle cell', 'fma226048'),
                                 ('Cytosol of smooth muscle cell', 'fma226046'),
                                 ('Cytosol of muscle cell', 'fma226044'), ('Cytosol of hemal cell', 'fma260695'),
                                 ('Cytosol of epithelial cell', 'fma260691')]

                    Description = []
                    for i, el in enumerate(valueRequired[singleSelection.value]):
                        if el == 'V_SS':
                            Description.append(widgets.HTML(value='Steady-state flux'))
                        elif el == 'V_m':
                            Description.append(widgets.HTML(value='Membrane voltage'))
                        else:
                            Description.append(widgets.HTML(value=[y[0] for y in entities if y[1] in
                                                                   [x for x in annots[singleSelection.value][el] if
                                                                    'entity$' in x][0]][0] + ' in ' +
                                                                  [y[0] for y in locations if y[1] in
                                                                   [x for x in annots[singleSelection.value][el] if
                                                                    'isPartOf$' in x][0]][0]))

                    # List of new headers
                    new_headers = copy.deepcopy(valueRequired[singleSelection.value])

                    # Dropdown widgets to select headers
                    dropdowns = []
                    selected_headers = {}

                    def create_dropdown(header):
                        dropdown = widgets.Dropdown(
                            options=[''] + df.columns.tolist(),
                            description=header + ':',
                            layout=Layout(width='250px', height='60px'),
                            style={'description_width': 'initial'},
                            disables=False)

                        # Apply the custom CSS class to the FloatText widget
                        dropdown.add_class('widget-dropdown')
                        dropdown.add_class('custom-text-style')
                        display(HTML(custom_css))

                        dropdowns.append(dropdown)

                        def dropdown_event(change, header=header):
                            selected_headers[header] = change.new

                        dropdown.observe(dropdown_event, names='value')

                    # Create dropdowns for each header
                    for header in new_headers:
                        create_dropdown(['$q_{}$'.format(header) if header not in ['V_SS', 'V_m'] else header][0] + [
                            '  ($u$mol/s)' if header == 'V_SS' else '  (mV)' if header == 'V_m' else '  ($u$mol)'][0])
                        display(HTML(custom_css))

                    # Button to show dropdowns
                    show_button = widgets.Button(description='Label Selection',
                                                 layout=Layout(width='300px', height='60px'))
                    show_button.style.font_family = 'monospace'
                    show_button.style.font_size = '28px'
                    show_button.style.button_color = 'lightblue'

                    dropdown_output = widgets.Output()

                    def show_button_clicked(button):
                        with dropdown_output:
                            dropdown_output.clear_output()
                            for i in range(len(Description)):
                                display(Description[i], dropdowns[i])
                                display(HTML(custom_css))
                            generate_button.layout.visibility = 'visible'

                    show_button.on_click(show_button_clicked)

                    # Button to generate new dataframe
                    generate_button = widgets.Button(description='Display the Steady-state Equation',
                                                     layout=Layout(width='600px', height='60px'))
                    generate_button.style.font_size = '28px'
                    generate_button.style.font_family = 'monospace'
                    generate_button.style.button_color = 'skyblue'
                    generate_button.layout.visibility = 'hidden'  # Initially hidden
                    generate_output = widgets.Output()

                    def generate_button_clicked(button):
                        with generate_output:
                            generate_output.clear_output()
                            selected_values = list(selected_headers.values())
                            if len(set(selected_values)) == len(new_headers) and all(selected_values):
                                if len(set(selected_values)) == len(selected_values):
                                    new_data = {['$q_{}$'.format(new_header) if new_header not in ['V_SS',
                                                                                                   'V_m'] else new_header][
                                                    0] + [
                                                    '  ($u$mol/s)' if new_header == 'V_SS' else '  (mV)' if new_header == 'V_m' else '  ($u$mol)'][
                                                    0]: df[selected_headers[['$q_{}$'.format(
                                        new_header) if new_header not in ['V_SS', 'V_m'] else new_header][0] + [
                                                                                '  ($u$mol/s)' if new_header == 'V_SS' else '  (mV)' if new_header == 'V_m' else '  ($u$mol)'][
                                                                                0]]] for new_header in new_headers}
                                    # global new_df
                                    new_df = pd.DataFrame(new_data)

                                    # global floatSSget
                                    floatSSget = []
                                    for header in new_df:
                                        floatSSget.append(widgets.FloatText(
                                            value=new_df[header].iloc[-1],
                                            description=header,
                                            disabled=False,
                                            layout=Layout(width='600px', height='60px')))

                                    for floatSS in floatSSget:
                                        if floatSS.description != 'V_m  (mV)':
                                            floatSS.value = copy.copy(floatSS.value / cellVol.value)

                                    for F in floatSSget:
                                        F.add_class('widget-dropdown')
                                        F.add_class('custom-text-style')

                                        # display(HTML(custom_css))

                                    pVals = pValsFunc(singleSelection, new_df, dataset=True)
                                    # Replace the parameter values in the equation
                                    equation = copy.deepcopy(equations[singleSelection.value])
                                    for param, value in pVals[singleSelection.value].items():
                                        param_pattern = re.escape('p_{' + param[1:] + '}')
                                        if abs(value) >= 0.01:
                                            formatted_value = "{:.2f}".format(value)
                                        else:
                                            formatted_value = "{:.2e}".format(value)
                                        equation = re.sub(param_pattern, formatted_value, equation)

                                    # Create a plot to display the expression
                                    fig, ax = plt.subplots(figsize=(6, 1))
                                    ax.text(1, 1, equation, fontsize=28, ha='center', va='center')
                                    ax.axis('off')
                                    plt.show()
                                else:
                                    print("Please select unique headers for all lists.")
                            else:
                                print("Please select headers for all lists (Check your selected data to be valid)")

                    generate_button.on_click(generate_button_clicked)

                    # Display buttons, dropdowns, and new dataframe output
                    display(show_button, dropdown_output, generate_button, generate_output)
                    display(HTML(custom_css))


                #################################################

                else:
                    print('Please upload a CSV file.')

        submit_button_widget.on_click(submit_button_clicked)
        vbox_result = widgets.VBox([file_upload_widget, submit_button_widget])
        display(vbox_result)


button_SS.on_click(button_SS_clicked)
buttonExperimentalData.on_click(buttonExperimentalData_clicked)

# Create the tab widgets
SS_page = widgets.VBox([button_SS, output_parameterization])
ED_page = widgets.VBox([buttonExperimentalData, output_csv])
display(HTML(custom_css))

# Create the Tab widget
tab = widgets.Tab()
tab.children = [SS_page, ED_page]
tab.set_title(0, "Option 1")
tab.set_title(1, "Option 2")


# Create a callback to clear the output when switching tabs
def clear_output_callback(change):
    output_parameterization.clear_output()
    output_csv.clear_output()


# Assign the callback to the selected_index attribute of the Tab widget
tab.observe(clear_output_callback, 'selected_index')

# Display the Tab widget and output
display(tab)
display(HTML(custom_css))

# # Section 4: Stimulus




button_stim = widgets.Button(
    description='Run',
    tooltip='Description',
    style={'description_width': 'initial'},
    button_style='warning',
    layout=Layout(width='120px', height='40px'))
button_stim.style.font_family = 'monospace'
button_stim.style.font_size = '28px'
outputStim = widgets.Output(layout={'border': '1px solid gray'})


def stimulate(event):
    with outputStim:
        outputStim.clear_output()
        img = Image.open('./figures/stimuli/step.jpg')
        display(img.resize((200, 200)))

        V_m = []
        if singleSelection.value == 'Template 5':
            V_m = [x.value for x in floatSSget if x.description == 'V_m' + '  (mV)']

        t_eval = np.arange(0, 100, 0.1)

        # Create radio buttons and slider
        global control_rb, control_slider, sliders

        control_rb = widgets.RadioButtons(
            options=[('q_{}'.format(x), x) for x in valueRequired1[singleSelection.value]] + [(x, x) for x in ['V_m'] if
                                                                                              V_m != []],
            description='Control variable:',
            style={'description_width': 'initial'},
            layout=Layout(width='100%'))

        control_slider = widgets.FloatSlider(
            value=[x.value for x in floatSSget if (x.description == '$q_{}$'.format(
                control_rb.value) + '  ($u$mol)' or x.description == 'V_m' + '  (mV)')][0], min=-200, max=1000,
            step=0.001,
            description='Control value:',
            style={'description_width': 'initial'},
            layout=Layout(width='100%'))

        sliders = {}
        for i, species in enumerate(valueRequired1[singleSelection.value] + [x for x in ['V_m'] if V_m != []]):
            sliders[species] = widgets.FloatSlider(
                value=[x.value for x in floatSSget if
                       x.description == '$q_{}$'.format(species) + '  ($u$mol)' or x.description == 'V_m  (mV)'][0],
                min=-200, max=1000, step=0.001,
                description=[x.description for x in floatSSget if
                             x.description == '$q_{}$'.format(species) + '  ($u$mol)' or x.description == 'V_m  (mV)'][
                                0] + ' :',
                layout=Layout(width='50%'))

        stepTime = FloatSlider(value=1000, min=0, max=1000, step=0.1, layout=Layout(width='100%'), readout_format='.3f',
                               style={'description_width': 'initial'}, description='t0')

        # Interact with the update function to update the plot
        widgets.interact(update_figure,
                         control_variable=control_rb,
                         control_value=control_slider,
                         stepTime=stepTime,
                         timespan=widgets.FloatSlider(value=1000, min=0, max=1000, step=0.1, description='Timespan:'),
                         **sliders
                         )

        # Save the selected conditions
        button_saveControl = widgets.Button(
            description='Save',
            tooltip='Save',
            style={'description_width': 'initial'},
            layout=Layout(width='120px', height='40px'))
        button_saveControl.style.font_family = 'monospace'
        button_saveControl.style.font_size = '20px'
        button_saveControl.style.button_color = 'lightgreen'
        output_saveControl = widgets.Output(layout={})

        def button_saveControl_clicked(event):
            with output_saveControl:
                clear_output()
                textControl = widgets.HTML(value="<h4><b>The selected control value has been saved.<b><h4>")
                display(textControl)
                for floatSS in floatSSget:
                    if '$q_{}$'.format(control_rb.value) + '  ($u$mol)' == floatSS.description:
                        floatSS.value = copy.copy(control_slider.value)

        button_saveControl.on_click(button_saveControl_clicked)
        vbox_result = widgets.VBox([button_saveControl, output_saveControl])
        pageSaveControl = widgets.HBox([vbox_result])
        display(pageSaveControl)


button_stim.on_click(stimulate)
vbox_result = widgets.VBox([button_stim, outputStim])
page7 = widgets.HBox([vbox_result])
display(page7)

# # Section 5: Save the Model




button_modelName = widgets.Textarea(
    description='Give it a name:',
    tooltip='Description',
    style={'description_width': 'initial'},
    layout=Layout(width='400px', height='50px'))

button_modelName.style.font_family = 'monospace'
button_modelName.style.font_size = '28px'
button_modelName.style.button_color = 'lightgreen'
button_modelName.add_class('widget-dropdown')
button_modelName.add_class('custom-text-style')
display(button_modelName)

button_addModel = widgets.Button(
    description='Save the model',
    tooltip='Description',
    style={'description_width': 'initial'},
    layout=Layout(width='250px', height='40px'))
button_addModel.style.font_family = 'monospace'
button_addModel.style.font_size = '28px'
button_addModel.style.button_color = 'lightgreen'
output9 = widgets.Output(layout={})


def button_addModel_clicked(event):
    with output9:
        output9.clear_output()
        if button_modelName.value.replace('\n', '') not in savedModels and button_modelName.value.replace('\n',
                                                                                                          '') != '' and button_modelName.value.replace(
                ' ', '') not in savedModels:
            savedModels.append(button_modelName.value.replace('\n', ''))
            text = widgets.HTML(value="<h3>Your generated model has been saved.</h3>")
            display(text)
            # Save the annotations as a new json file
            f = open("./temporary files/selectedAnnotations.json")
            file = json.load(f)
            file[button_modelName.value.replace('\n', '')] = file[singleSelection.value]
            del file[singleSelection.value]
            savedAnnotations.update(file)
            with open('./savedAnnotations.json', 'w') as f:
                json.dump(savedAnnotations, f, indent=4, sort_keys=True, ensure_ascii=False)

                # Save the P values as a new json file
            f = open("./temporary files/pVals.json")
            file = json.load(f)
            file[button_modelName.value.replace('\n', '')] = file[singleSelection.value]
            del file[singleSelection.value]
            savedpVals.update(file)
            with open('./savedModels/savedpVals.json', 'w') as f:
                json.dump(savedpVals, f, indent=4, sort_keys=True, ensure_ascii=False)

                # Save the type as a new json file
            file = dict([(button_modelName.value.replace('\n', ''), singleSelection.value)])
            savedtypes.update(file)
            with open('./savedModels/savedtypes.json', 'w') as f:
                json.dump(savedtypes, f, indent=4, sort_keys=True, ensure_ascii=False)


        elif button_modelName.value.replace('\n', '') == '':

            text = widgets.HTML(value="<h3><font color='red'>None: Enter a valid name.</h3>")
            display(text)
        else:
            text = widgets.HTML(value="<h3><font color='red'>Repetitive name: Enter a new name.</h3>")
            display(text)

        print('Your saved model is:')
        modelUnits = {'Concentration': 'micro', 'Flow rate': 'milli'}
        model = singleModelCellmlGenerator(modelUnits, savedpVals, button_modelName.value)
        xmlAnnot(model, savedAnnotations)


button_addModel.on_click(button_addModel_clicked)
vbox_result = widgets.VBox([button_addModel, output9])
display(vbox_result)

# <div style="text-align: center;">
#     <h3>👉🏻 Upload the saved model on a PMR workspace 👈🏻</h3>
#     <img src="./figures/physiome.png" alt="Figure" width="200" height="200">
# </div>

# <div class="alert alert-block alert-warning">
#     <h1><b>Model Composition</b></h1>
# </div>
#
# ***

# # Section 1: Search on PMR




search_modelName = widgets.Textarea(
    description='Enter a model name:',
    tooltip='Description',
    style={'description_width': 'initial'},
    layout=Layout(width='450px', height='50px')
)
search_modelName.style.font_family = 'monospace'
search_modelName.style.font_size = '28px'
search_modelName.style.button_color = 'red'
search_modelName.add_class('widget-dropdown')
search_modelName.add_class('custom-text-style')
display(search_modelName)

button_PMR = widgets.Button(
    description='Search on PMR',
    tooltip='Search',
    style={'description_width': 'initial'},
    button_style='info',
    layout=Layout(width='300px', height='40px'))

button_PMR.style.font_family = 'monospace'
button_PMR.style.font_size = '28px'
button_PMR.style.button_color = 'red'
button_PMR.add_class('widget-dropdown')
button_PMR.add_class('custom-text-style')
outputPMR = widgets.Output()

global downloadedModels, xmlParticles, xmlGrouped
downloadedModels = []
xmlParticles = {};
xmlGrouped = {}

dir_name = "./savedModels/PMR"
test = os.listdir(dir_name)

for item in test:
    if item.endswith(".cellml"):
        os.remove(os.path.join(dir_name, item))


def button_PMR_clicked(event):
    with outputPMR:
        clear_output()

        pmrModel = pmrSearching(search_modelName.value.replace(' ', '%20'))

        xmlParticles = RDFpmrSearching(search_modelName.value.replace(' ', '%20'))

        buttons = []
        links = []
        outputs = []

        if pmrModel != []:
            for i, el in enumerate(pmrModel):
                links.append(el)

            for link in links:
                text = widgets.HTML(value=link)
                button = widgets.Button(description='Select', tooltip='Select', style={'description_width': 'initial'},
                                        button_style='info')
                output = widgets.Output()
                display(widgets.HBox([text, button]))
                display(output)
                buttons.append(button)
                outputs.append(output)

            def on_button_click(button):
                with outputs[buttons.index(button)]:
                    clear_output()
                    text = widgets.HTML("<h3><font color='red'>Selected!</h3>")
                    display(text)
                    filename = os.path.basename(links[buttons.index(button)]).split('.')[0]
                    folder_path = './savedModels/PMR'
                    # Check if the folder exists, and create it if it doesn't
                    if not os.path.exists(folder_path):
                        os.makedirs(folder_path)

                    file = "./savedModels/PMR/{}.cellml".format(filename)
                    urllib.request.urlretrieve(links[buttons.index(button)], file)
                    downloadedModels.append(filename + '.cellml')

                    for workspace in xmlParticles:
                        for ID, predicate, obj in xmlParticles[workspace]:
                            if workspace in links[buttons.index(button)] and ID.split('#')[0] in links[
                                buttons.index(button)]:
                                if ID.split('#')[0] not in xmlGrouped:
                                    xmlGrouped[ID.split('#')[0]] = {}  # model name
                                    if ID.split('#')[1] not in xmlGrouped[ID.split('#')[0]]:
                                        xmlGrouped[ID.split('#')[0]][ID.split('#')[1]] = []
                                        xmlGrouped[ID.split('#')[0]][ID.split('#')[1]].append((predicate, obj))
                                    else:
                                        xmlGrouped[ID.split('#')[0]][ID.split('#')[1]].append((predicate, obj))
                                else:
                                    if ID.split('#')[1] not in xmlGrouped[ID.split('#')[0]]:
                                        xmlGrouped[ID.split('#')[0]][ID.split('#')[1]] = []
                                        xmlGrouped[ID.split('#')[0]][ID.split('#')[1]].append((predicate, obj))
                                    else:

                                        xmlGrouped[ID.split('#')[0]][ID.split('#')[1]].append((predicate, obj))

            for button in buttons:
                button.on_click(lambda b: on_button_click(b))


        else:
            print('No such a file on PMR')


button_PMR.on_click(button_PMR_clicked)
vbox_result = widgets.VBox([button_PMR, outputPMR])
display(vbox_result)

# # Section 2: Model Composition




button_compose = widgets.Button(
    description='Compose the models',
    tooltip='Description',
    style={'description_width': 'initial'},
    layout=Layout(width='300px', height='40px'))
button_compose.style.font_family = 'monospace'
button_compose.style.font_size = '28px'
button_compose.style.button_color = 'orange'
compositionOutput = widgets.Output(layout={})


def button_compose_clicked(event):
    with compositionOutput:
        compositionOutput.clear_output()

        PMRmodelComposition()

        file_name = 'CompositeModel'
        pythonCreator('./', file_name)

        # get_ipython().run_line_magic('run', '{file_name}.py')
        subprocess.run(["python", {file_name}.py])

        ##############################################################################################################################################
        my_rates = np.zeros((STATE_COUNT))
        my_variables = np.zeros((VARIABLE_COUNT))
        my_state_variables = np.zeros((STATE_COUNT))

        initialise_variables(my_state_variables, my_variables)

        # Create sliders
        global control_slider
        control_slider = {}
        for i, species in enumerate(speciesNoDuplicate):
            control_slider[species[-1]] = widgets.FloatSlider(
                value=my_state_variables[i], min=0,
                max=[my_state_variables[i] * 10 if my_state_variables[i] != 0 else 10][0], step=0.001,
                description=species[-1] + ' :',
                layout=Layout(width='50%'))

        stepSize = widgets.Dropdown(options=[0.001, 0.01, 0.1, 1.0, 10, 100],
                                    value=0.001,
                                    description='Step size:',
                                    style={'description_width': 'initial'},
                                    disabled=False)

        # timespan=widgets.FloatSlider(value=10, min=0, max=100000, step=stepSize.value, description='Timespan (ms) :', style={'description_width': 'initial'}, layout=Layout(width='50%'))
        timespan_text = widgets.FloatText(value=10, style={'description_width': 'initial'},
                                          description='Timespan (ms):')

        # Define the update function
        def update_composed_figure(T, stepSize, **control_value):

            my_rates = np.zeros((STATE_COUNT))
            my_variables = np.zeros((VARIABLE_COUNT))
            my_state_variables = np.zeros((STATE_COUNT))

            initialise_variables(my_state_variables, my_variables)
            for i, species in enumerate(speciesNoDuplicate):
                my_state_variables[i] = control_value[species[-1]]
            compute_computed_constants(my_variables)
            compute_rates('t', my_state_variables, my_rates, my_variables)

            global Solution
            step_size = stepSize
            step_count = int(T / step_size)
            Solution = {}

            for step in range(0, step_count):
                time = step * step_size

                Solution[step] = []
                Solution[step].append(time)

                # Compute the rates at this step using the given function.
                compute_rates(time, my_state_variables, my_rates, my_variables)

                # Compute the states.
                for s in range(0, STATE_COUNT):
                    my_state_variables[s] = my_state_variables[s] + my_rates[s] * step_size

                # Compute the variables.
                compute_variables(time, my_state_variables, my_rates, my_variables)

                for s in range(0, VARIABLE_COUNT):
                    Solution[step].append(my_variables[s])

                for s in range(0, STATE_COUNT):
                    Solution[step].append(my_state_variables[s])

            # Update the plot
            entities = [('H+', 'CHEBI:15378'), ('HCO3-', 'CHEBI:17544'), ('K+', 'CHEBI:29103'), ('Na+', 'CHEBI:29101'),
                        ('Mg2+', 'CHEBI:18420'), ('Cl-', 'CHEBI:17996'), ('Ca2+', 'CHEBI:29108'),
                        ('Fe2+', 'CHEBI:29033'), ('P', 'CHEBI:30207')]
            locations = [('Extracellular environment', 'fma70022'), ('Cytosol of stem cell', 'fma260697'),
                         ('Cytosol of neuron', 'fma226054'), ('Cytosol of connective tissue cell', 'fma260699'),
                         ('Cytosol of skeletal muscle fiber', 'fma226052'), ('Cytosol of cardiac myocyte', 'fma226050'),
                         ('Cytosol of striated muscle cell', 'fma226048'),
                         ('Cytosol of smooth muscle cell', 'fma226046'), ('Cytosol of muscle cell', 'fma226044'),
                         ('Cytosol of hemal cell', 'fma260695'), ('Cytosol of epithelial cell', 'fma260691')]

            matplotlib.rcParams.update({'font.size': 15})
            fig, ax = plt.subplots(nrows=[int(len(speciesNoDuplicate) / 2) if len(speciesNoDuplicate) % 2 == 0 else int(
                len(speciesNoDuplicate) / 2) + 1][0], ncols=2, figsize=(25, 7 * len(speciesNoDuplicate) / 2))
            plt.subplots_adjust(
                wspace=0.2,
                hspace=0.4)
            for i, species in enumerate(speciesNoDuplicate):
                line = []
                timecourse = []

                for step in Solution:
                    timecourse.append(list(Solution[step])[0])
                    line.append(list(Solution[step])[-(len(speciesNoDuplicate) - i)])

                if i < [int(len(speciesNoDuplicate) / 2) if len(speciesNoDuplicate) % 2 == 0 else int(
                        len(speciesNoDuplicate) / 2) + 1][0]:
                    k = copy.copy(i)
                    j = 0
                else:
                    k = copy.copy(i - [int(len(speciesNoDuplicate) / 2) if len(speciesNoDuplicate) % 2 == 0 else int(
                        len(speciesNoDuplicate) / 2) + 1][0])
                    j = 1

                ax[k, j].plot(timecourse, line, color=matplotlib.cm.hsv(float(i) / len(speciesNoDuplicate)),
                              linewidth=5, label=species[-1])
                ax[k, j].set_xlabel('Time (s)')
                ax[k, j].grid(color='gray', linestyle='-', linewidth=0.5)
                ax[k, j].set_ylabel('Amount ($u$mol)')
                ax[k, j].legend(loc='best')
                ax[k, j].patch.set_facecolor('black')

                for x in species:
                    for y in locations:
                        if x[0] == 'http://biomodels.net/biology-qualifiers/isPartOf' and y[1] in x[1]:
                            location = [y[0]]
                    for y in entities:
                        if y[1] in x[1] and x[
                            0] == 'http://biomodels.net/biology-qualifiers/isVersionOf' and 'https://identifiers.org/opb' not in \
                                x[1]:
                            entity = [y[0]]

                if entity != [] and location != []:
                    title = entity + [' in '] + location
                if entity == [] and location != []:
                    title = ['unidentified entity in '] + location
                if location == [] and entity != []:
                    title = entity + [' in unidentified location']
                if entity == [] and location == []:
                    title = ['unidentified entity'] + [' in unidentified location']
                ax[k, j].set_title(''.join(title))

            plt.show()

        # Interact with the update function to update the plot
        widgets.interact(update_composed_figure,
                         T=timespan_text,
                         stepSize=stepSize,
                         **control_slider
                         )

        button_saveSS = widgets.Button(
            description='Save',
            tooltip='Save',
            style={'description_width': 'initial'},
            layout=Layout(width='120px', height='40px'))
        button_saveSS.style.font_family = 'monospace'
        button_saveSS.style.font_size = '20px'
        button_saveSS.style.button_color = 'lightgreen'
        outputInteract = widgets.Output(layout={})

        def button_saveSS_clicked(event):
            with outputInteract:
                clear_output()
                text = widgets.HTML(value="<h3><b>The model is updated with the selected values.<b><h3>")
                display(text)
                modelUpdate = parse_model(os.path.join('./', 'CompositeModel.cellml'), False)
                for key in control_slider.keys():
                    for compNum in range(modelUpdate.componentCount()):
                        for varNum in range(modelUpdate.component(compNum).variableCount()):
                            if modelUpdate.component(compNum).variable(varNum).name() == key:
                                modelUpdate.component(compNum).variable(varNum).setInitialValue(
                                    control_slider[key].value)
                printer = libcellml.Printer()
                # print(printer.printModel(modelUpdate))
                writeCellML(modelUpdate, printer, 'CompositeModel')

        button_saveSS.on_click(button_saveSS_clicked)

        vbox_result = widgets.VBox([button_saveSS, outputInteract])
        pageInteract = widgets.HBox([vbox_result])
        display(pageInteract)


##########################################################################################################################################################


button_compose.on_click(button_compose_clicked)
vbox_result = widgets.VBox([button_compose, compositionOutput])
display(vbox_result)

# # Section 3: Graph Representation




button_graph = widgets.Button(
    description='Plot the graph',
    tooltip='Description',
    style={'description_width': 'initial'},
    layout=Layout(width='300px', height='40px'))
button_graph.style.font_family = 'monospace'
button_graph.style.font_size = '28px'
button_graph.style.button_color = 'gray'
button_graph.add_class('widget-dropdown')
button_graph.add_class('custom-text-style')
graphOutput = widgets.Output(layout={})

figSize = widgets.FloatText(value=20, description='Figure size: ', layout=Layout(width='250px', height='40px'))
figSize.add_class('widget-dropdown')
figSize.add_class('custom-text-style')

display(figSize)


def button_graph_clicked(event):
    with graphOutput:
        graphOutput.clear_output()

        graphPlot(newCellmlNames, Species, speciesNoDuplicate, figSize.value)


button_graph.on_click(button_graph_clicked)
vbox_result = widgets.VBox([button_graph, graphOutput])
display(vbox_result)








