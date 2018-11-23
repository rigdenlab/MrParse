'''
Created on 23 Nov 2018

@author: jmht
'''


NULL_SYMBOL = '-'

class AnnotationSymbol(object):
    def __init__(self):
        self.symbol = None
        self.name = None
        self.colour = None
        self.parent = None

    
class SequenceAnnotation(object):
    def __init__(self):
        self.source = None
        self.probabilties = None
        self.annotation = None
        self.symbols = dict()
        
    def add_symbol(self, sobj):
        self.symbols[sobj.symbol] = sobj
        sobj.parent = self
    
    @property
    def length(self):
        return len(self.annotation)

    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str
        

class AnnotationChunk(object):
    def __init__(self, start=None, end=None, annotation=None):
        self.start = start
        self.end = end
        self.annotation = annotation
    
    @property
    def stype(self):
        return self.annotation.symbol
    
    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str
    
    
def get_annotation_chunks(annotation, annotation_data):
    chunks = []
    chunk = None
    symbol_map = {}
    for a in annotation_data:
        symbol_map.update(a.symbols)
    for i, s in enumerate(annotation):
        if s in symbol_map.keys():
            if not chunk:
                chunk = AnnotationChunk(start=i, annotation=symbol_map[s])
            elif chunk.stype != s:
                chunk.end = i
                chunks.append(chunk)
                chunk = AnnotationChunk(start=i, annotation=symbol_map[s])
        else:
            if chunk:
                chunk.end = i
                chunks.append(chunk)
                chunk = None
    if chunk:
        chunk.end = i
        chunks.append(chunk)
    return chunks