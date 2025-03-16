import logging

class Chi2Ext:
    def __init__(self):
        self.logger = logging.getLogger(f'{__name__}.{__class__.__name__}')
        self.logger.debug('Hooray, instantiated...')

    def get_chi2(self, parset):
        self.logger.debug(f'This is get_chi2 and I got {parset=}.')
        return 27
