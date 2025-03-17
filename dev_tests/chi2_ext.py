import logging

class Chi2Ext:
    def __init__(self, **kwds):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')

        args = tuple(f'{a}={kwds[a]}' for a in kwds)
        self.logger.debug(f'Instantiated with parameters {args}.')

    def chi2(self, parset):
        self.logger.debug(f'This is get_chi2 and I got {parset=}.')
        return 27
