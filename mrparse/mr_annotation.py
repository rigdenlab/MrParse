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

    
class SequenceAnnotation(object):
    def __init__(self):
        self.source = None
        self.probabilties = None
        self.annotation = None
        self.annotation_symbols = None
    
    @property
    def length(self):
        return len(self.annotation)

    @property
    def symbols(self):
        return [s.symbol for s in self.annotation_symbols]

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
    
    def __str__(self):
        attrs = [k for k in self.__dict__.keys() if not k.startswith('_')]
        INDENT = "  "
        out_str = "Class: {}\nData:\n".format(self.__class__)
        for a in sorted(attrs):
            out_str += INDENT + "{} : {}\n".format(a, self.__dict__[a])
        return out_str
    
    
def get_annotation_chunks(annotation):
    chunks = []
    chunk = None
    _annotations = { a.symbol : a for a in annotation }
    for i, s in enumerate(annotation.annotation):
        if s in _annotations.keys():
            if not chunk:
                chunk = AnnotationChunk(start=i, annotation=_annotations[s])
            elif chunk.stype != s:
                chunk.end = i
                chunks.append(chunk)
                chunk = AnnotationChunk(start=i, annotation=_annotations[s])
        else:
            if chunk:
                chunk.end = i
                chunks.append(chunk)
                chunk = None
    if chunk:
        chunk.end = i
        chunks.append(chunk)
    return chunks