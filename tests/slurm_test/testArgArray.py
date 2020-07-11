#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 13:15:33 2020

@author: twguest
"""


import sys

param = sys.argv[1]

def main(variable):
    
    var = int(variable)*2
    
    return var

if __name__ == '__main__':
    var = main(param)
    print(var)    
    
    