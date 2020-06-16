class Result(object):
    '''
    Convert dictionary to a class, allowing easier access, i.e.
        dict.item    vs.    dict['item']
    '''
    def __init__(self, d):
        self.__dict__ = d
