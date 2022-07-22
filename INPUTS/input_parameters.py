#!/usr/bin/env python3

import yaml

qgtMin = {'qgtMin':0.00, 'Units':'m^3/s', 'Definition': 'minimum inlet gas flow rate'}

qgtMax = {'qgtMax':0.00, 'Units':'m^3/s', 'Definition': 'maximum inlet gas flow rate'}

qgtIncrement = {'qgtIncrement':0.00, 'Units':'m^3/s', 'Definition': 'step size for inlet gas flow rate'}

qlvr = {'qlvr': 0.00, 'Units': 'm^3/s', 'Definition': 'inlet liquid flow rate'}

rhop = {'rhop': 0.00, 'Units': 'kg/m^3', 'Definition': 'particle density'}

rhol = {'rhol': 0.00, 'Units': 'kg/m^3', 'Definition': 'liquid density'}

rhog = {'rhog': 0.00, 'Units': 'kg/m^3', 'Definition': 'gas density'}

mul = {'mul': 0.00, 'Units': 'Pa.s', 'Definition': 'liquid viscosity'}

sigma = {'sigma': 0.00, 'Units': 'N/m', 'Definition': 'gas-liquid interfacial tension'}

P = {'P': 0.00, 'Units': 'MPa', 'Definition': 'operating pressure'}

mcat = {'mcat': 0.00, 'Units': 'kg', 'Definition': 'mass of catalyst inventory'}

Vp = {'Vp': 0.00, 'Units': 'm^3', 'Definition': 'particle volume'}

Ap = {'Ap': 0.00, 'Units': 'm^2', 'Definition': 'particle surface area'}

Dc = {'Dc': 0.00, 'Units': 'm', 'Definition': 'column diameter'}

drec = {'drec': 0.00, 'Units': 'm', 'Definition': 'recycle line diameter'}

vsep = {'vsep': 0.00, 'Units': 'm^3', 'Definition': 'gas-liquid separator volume'}

dor = {'dor': 0.00, 'Units': 'm', 'Definition': 'distributor (i.e., bubble cap grid) outlet orifice diameter'}

nor = {'nor': 0.00, 'Units': 'non-dimensional', 'Definition': 'number of outlet orifices per riser'}

nr = {'nr': 0.00, 'Units': 'non-dimensional', 'Definition': 'number of risers'}

classWidth = {'classWidth': 0.00, 'Units': 'mm', 'Definition': 'bubble class width'}

minCutOff = {'minCutOff': 0.00, 'Units': 'non-dimensional', 'Definition': 'minimum cumulative probability (for constructing bubble classes)'}

maxCutOff = {'maxCutOff': 0.00, 'Units': 'non-dimensional', 'Definition': 'maximum cumulative probability (for constructing bubble classes)'}

beta1 = {'beta1': 0.00, 'Units': '1/s', 'Definition': 'geometry dependent parameter of gas-liquid separation submodel'}

beta2 = {'beta2': 0.00, 'Units': '1/mm', 'Definition': 'bubble and geometry dependent parameter of gas-liquid separation submodel'}

Rstep = {'Rstep': 0.00, 'Units': 'non-dimensional', 'Definition': 'step size for the liquid recycle ratio'}

hset = {'hset': 0.00, 'Units': 'm', 'Definition': 'desired bed height'}

absTol = {'absTol': 0.00, 'Units': 'problem dependent', 'Definition': 'absolute tolerance'}

relTol = {'relTol': 0.00, 'Units': 'non-dimensional', 'Definition': 'relative tolerance'}

maxIter = {'maxIter': 0, 'Units': 'non-dimensional', 'Definition': 'maximum number of iterations'}

parameters = [qgtMin, qgtMax, qgtIncrement, qlvr, rhop, rhol, rhog, mul, sigma, P, mcat, Vp, Ap, Dc,
              drec, vsep, dor, nor, nr, classWidth, minCutOff, maxCutOff, beta1, beta2, Rstep, hset, 
              absTol, relTol, maxIter]

with open('input_parameters.yaml', 'w') as f:
    data = yaml.dump(parameters, f, sort_keys=False)
