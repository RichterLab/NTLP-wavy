#!/usr/bin/env python3
import os
import sys

# Copies all files from first arg to second arg (dir created if needed)
# Changes param directories to /scratch365/$USER/[this dir, 2 deep]/*
# Also sets les.run process name to the current directory, 2 deep


def setPath(text, var, path):
    home = os.path.expanduser("~")
    user = home.split('/')[-1]
    dir = "/scratch365/" + user + "/"
    i = text.find("\n&path_names")
    i = text.find(var, i)
    i = text.find('\"', i) + 1
    j = text.find('\"', i)
    j = text[i:j].rfind('/') + i
    return text[:i] + dir+path + text[j:]

def setProcName(text, name):
    var = "#$ -N "
    i = text.find(var)+len(var)
    j = text.find('\n', i)
    return text[0:i] + name + text[j:]

def setVar(text, var, val):
    # Find variable
    i = text.find(var) + len(var);
    if i == len(var)-1: return text # var not found
    # Skip '=' and spaces
    while text[i]==' ' or text[i]=='=': i+=1
    # Find end of current value
    j = i
    while text[j] != '!' and text[j] != '\n': j+=1
    # Insert new value and return
    return text[:i] + str(val) + text[j:]

def setArg(text, prog, arg, val):
    # Find command
    i = text.find(prog) + len(prog)
    # Skip to arg
    i = text.find(arg, i) + len(arg)
    # Skip spaces
    while text[i]==' ': i+=1
    j = i
    # Find end of current value
    while text[j]!=" " and text[j]!=":" and text[j]!="\n": j+=1
    # Insert new value and return
    return text[:i] + str(val) + text[j:]

def printHelp():
    print("Usage:\n  ./"+sys.argv[0][1:-3]+" oldRun newRun -param1 value -param2 value " \
        + "-paramN value")
    print("Help:\n  Copies oldRun directory to newRun directory, which is created. \n" \
        + "  The paths set in params.in are updated to /scratch365/$USER/[LESdir]/newRun \n" \
        + "  and the process name is set to newRun_[LESdir] in les.run. Optional \n" \
        + "  parameters matching entries in params.in are set to the given values")



# Parse CLI #
args = {} # dictionary to hold cli args
for i in range(0, len(sys.argv)):
    a = sys.argv[i]
    if a == '-help' or a =='-h':
        printHelp()
        quit()
    elif a[0] == '-':
        args[a[1:]] = sys.argv[i+1]

if len(sys.argv) < 3:
    print("Error:\n  Not enough arguments")
    printHelp()
    quit()

# Open files #
oldRun = sys.argv[1]
newRun = sys.argv[2]
with open(oldRun+"/params.in") as file:
    pStr = file.read()
with open(oldRun+"/les.run") as file:
    lStr = file.read()


# Set parameters given
for a in args:
    pStr = setVar(pStr, a, args[a] + "   ")

# Get parent/newRun dir for saving to same dir in scratch
dirs = os.getcwd().split('/')
procName = newRun + "_" + dirs[-1]
path = dirs[-1] + "/" + newRun

# Copy over run files
os.system("mkdir " + newRun + " &> /dev/null")  # silences output
os.system("cp " + oldRun+"/*.dat " + newRun+"/") # all dat files
os.system("cp " + oldRun+"/clean " + newRun+"/clean")
os.system("cp " + oldRun+"/report " + newRun+"/report")

# Change path in params.in
pStr = setPath(pStr, "path_seed", path)
pStr = setPath(pStr, "path_part", path)
pStr = setPath(pStr, "path_res",  path)
pStr = setPath(pStr, "path_ran",  path)

# Change process name in les.run
lStr = setProcName(lStr, procName)

with open(newRun + "/params.in", "w+") as file:
    file.write(pStr)
with open(newRun + "/les.run", "w+") as file:
    file.write(lStr)
