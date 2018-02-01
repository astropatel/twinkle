
import os
homedir = os.path.expanduser('~')

def getDir(dir):
    for root, dirs, files in os.walk(homedir):
        if dir in dirs:
            return os.path.join(root, dir)

def SupportFiles():
    return getDir('SupportFiles')

def NextGen():
    return getDir('NextGen')

def Atlas9():
    return getDir('Atlas9')

def RSR():
    return getDir('RSR')

def WorkingDir(name):
    return getDir(name)