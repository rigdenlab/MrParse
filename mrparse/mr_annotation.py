'''
Created on 23 Nov 2018

@author: jmht
'''
import copy


class AnnotationSymbol(object):
    def __init__(self, name=None, symbol=None, stype=None):
        __slots__ = ('name', 'symbol', 'stype', 'score', 'source')
        self.name = name
        self.symbol = symbol
        self.stype = stype
        self.score = None
        self.source = None

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


NULL_ANNOTATION = AnnotationSymbol()
NULL_ANNOTATION.name ='null'
NULL_ANNOTATION.symbol ='-'
NULL_ANNOTATION.stype = 'null'
NULL_ANNOTATION.score = 0.0
NULL_ANNOTATION.source = 'null'



class SequenceAnnotation(object):
    def __init__(self, null_symbol=NULL_ANNOTATION.symbol):
        __slots__ = ('source', 'scores', 'annotation', 'annotation_library', 'null_symbol')
        self.source = None
        self.scores = []
        self.annotation = '' # list of annotation symbols
        self.annotation_library = dict()
        self.null_symbol = null_symbol
    
    def add_annotation(self, annotation):
        if annotation != NULL_ANNOTATION:
            assert self.annotation_is_significant(annotation), "Cannot find: %s" % annotation
        self.annotation += annotation.symbol
        self.scores.append(annotation.score)
    
    def annotation_is_significant(self, annotation):
        return annotation is not NULL_ANNOTATION and annotation in self.annotation_library.values()
    
    def has_annotation_symbol(self, symbol):
        return symbol in self.annotation_library
    
    def library_add_annotation(self, annotation):
        assert isinstance(annotation, AnnotationSymbol)
        annotation.source = self.source
        self.annotation_library[annotation.symbol] = annotation
        
    def __getitem__(self, idx):
        symbol = self.annotation[idx]
        if self.has_annotation_symbol(symbol):
            a = copy.copy(self.annotation_library[symbol])
            a.score = self.scores[idx]
#         elif symbol == self.null_symbol:
        else:
            a = copy.copy(NULL_ANNOTATION)
            a.symbol = self.null_symbol
        return a

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError(other)
        assert len(self) == len(other)
        ca = self.__class__() # Class holding consensus annotation
        
        # Combine annotation_libraries - is it worth removing those entries that aren't
        # included in the consensus?
        ca.annotation_library = dict(self.annotation_library, **other.annotation_library)
        
        # For now just use prediction and leave probabiltiies
        for i in range(len(self)):
            if self.annotation_is_significant(self[i]) and other.annotation_is_significant(other[i]):
                ca.add_annotation(NULL_ANNOTATION)
            elif self.annotation_is_significant(self[i]):
                ca.add_annotation(self[i])
            elif other.annotation_is_significant(other[i]):
                ca.add_annotation(other[i])
            else:
                ca.add_annotation(NULL_ANNOTATION)
        return ca        

    def __len__(self):
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
    for i, a in enumerate(annotation):
        if a != NULL_ANNOTATION:
            if not chunk:
                chunk = AnnotationChunk(start=i, annotation=a)
            elif chunk.annotation != a:
                chunk.end = i
                chunks.append(chunk)
                chunk = AnnotationChunk(start=i, annotation=a)
        else:
            if chunk:
                chunk.end = i
                chunks.append(chunk)
                chunk = None
    if chunk:
        chunk.end = i
        chunks.append(chunk)
    return chunks