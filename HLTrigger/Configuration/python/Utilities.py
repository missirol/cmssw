import HLTrigger.Configuration.Tools.options as _options
import HLTrigger.Configuration.Tools.confdb as _confdb

def _build_options(**args):
    options = _options.HLTProcessOptions()
    for key, val in args.items():
        setattr(options, key, val)

    return options


def getHltConfiguration(menu, **args):
    args['menu'] = menu
    args['fragment'] = False
    options = _build_options(**args)

    hlt = {'process': None, 'fragment': None}
    exec(_confdb.HLTProcess(options).dump(), globals(), hlt)
    hlt = hlt['process'] if hlt['process'] != None else hlt['fragment']

    return hlt


def loadHltConfiguration(process, menu, **args):
    args['menu'] = menu
    args['fragment'] = True
    options = _build_options(**args)

    hlt = {'process': None, 'fragment': None}
    exec(_confdb.HLTProcess(options).dump(), globals(), hlt)
    hlt = hlt['process'] if hlt['process'] != None else hlt['fragment']

    process.extend( hlt )

    return process


import FWCore.ParameterSet.Config as _cms
_cms.Process.loadHltConfiguration = loadHltConfiguration
