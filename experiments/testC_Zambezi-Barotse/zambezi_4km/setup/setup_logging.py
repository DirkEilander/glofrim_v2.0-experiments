# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 09:38:02 2019

@author: Arjen Haag (haag)
"""

import sys
import logging


# main logging setup
def setupLogging(logPath, name='default'):
    """
    Setup/config of logging procedures.
    
    Parameters:
        logPath : [string] path specifying where to store logging file
            
    Output: instance of Python logging module extended with counting functionality
    """
    # setup basic logging configuration
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)-8.8s]  %(message)s",
        handlers=[
            logging.FileHandler("{0}/{1}.log".format(logPath, 'setup')),
            logging.StreamHandler()
        ])
    # decorator to determine number of calls for a method (https://stackoverflow.com/questions/812477/how-many-times-was-logging-error-called)
    class CallCounted:
        def __init__(self,method):
            self.method=method
            self.counter=0
        def __call__(self,*args,**kwargs):
            self.counter+=1
            return self.method(*args,**kwargs)
    # create a logger
    logger = logging.getLogger(name)
    # extend with decorator
    logger.info    = CallCounted(logger.info)
    logger.warning = CallCounted(logger.warning)
    logger.error   = CallCounted(logger.error)
    logger.fatal   = CallCounted(logger.fatal)
    return logger


# close logger (so a clean one can be started for a new catchment)
def closeLogging():
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)


# functions for logging at specific levels, or printing to console if no logger was configured
def showInfo(logger, mssg):
    if logger != None:
        logger.info(mssg)
    else:
        print(mssg)
def showWarning(logger, mssg):
    if logger != None:
        logger.warning(mssg)
    else:
        print('WARNING! ' + mssg)
def showError(logger, mssg):
    if logger != None:
        logger.error(mssg)
    else:
        print('ERROR! ' + mssg)
def showFatal(logger, mssg, exit=True):
    if logger != None:
        logger.fatal(mssg)
    else:
        print('ERROR! ' + mssg)
    if exit:
        sys.exit()


# statistics
def getStats(logger):
    log_warns  = logger.warning.counter
    log_errors = logger.error.counter
    log_fatal = logger.fatal.counter
    logger.info('Checking logging statistics...')
    logger.info('warning:  ' + str(log_warns))
    logger.info('error:    ' + str(log_errors))
    logger.info('fatal:    ' + str(log_fatal))
    return {'warn':log_warns, 'err':log_errors, 'crit':log_fatal}