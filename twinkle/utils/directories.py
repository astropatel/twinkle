
__author__ = 'Joe Trollo'

import os
homedir = os.path.expanduser('~')

def getDir(dir):
    for root, dirs, files in os.walk(homedir):
        if dir in dirs:
            return os.path.join(root, dir)

def SupportFiles():
    """

    Returns:

    """
    return getDir('Inputs_and_Models')

def NextGen():
    """

    Returns

    """
    return getDir('NextGen')

def Atlas9():
    """

    Returns:


    """
    return getDir('Atlas9')

def RSR():
    """

    Returns:

    """
    return getDir('RSR')

def WorkingDir(name):
    """

    Args:
        name

    Returns:


    """
    return getDir(name)