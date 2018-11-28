'''
Created on 23 Nov 2018

@author: jmht
'''
import copy


class AnnotationSymbol(object):
    def __init__(self, name=None, symbol=None, stype=None):
        __slots__ = ('name', 'symbol', 'stype', 'score')
        self.name = name
        self.symbol = symbol
        self.stype = stype
        self.score = None


NULL_ANNOTATION = AnnotationSymbol()
NULL_ANNOTATION.name ='null'
NULL_ANNOTATION.symbol ='-'
NULL_ANNOTATION.stype ='null'


class SequenceAnnotation(object):
    def __init__(self):
        __slots__ = ('source', 'scores', 'annotations', 'annotation_library')
        self.source = None
        self.scores = None
        self.annotations = None
        self.annotation_library = dict()
    
    def library_add_annotation(self, annotation):
        assert isinstance(annotation, AnnotationSymbol)
        self.annotation_library[annotation.symbol] = annotation
        
    def __getitem__(self, idx):
#         if not isinstance(idx, int) :
#             raise TypeError(idx)
#         if self.annotations is not None and 0 < idx < len(self.annotations):
#             raise IndexError(idx)
#         assert len(self.annotations) == len(self.probabilties)
        symbol = self.annotations[idx]
        a = copy.copy(self.annotation_library[symbol])
        a.score = self.scores[idx]
        return a

    def __len__(self):
        return len(self.annotations)
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
#             return self.__dict__ == other.__dict__
            return other.stype == self.stype
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

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