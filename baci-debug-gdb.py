#
# STK pretty printers for use with gdb 7.0.1 on x86 linux, gcc compiled 
#
# see http://tromey.com/blog/?p=494 for explanations
#

import gdb
import sys
sys.path.append("/usr/share/gcc-4.5/python/libstdcxx")

from v6.printers import register_libstdcxx_printers
register_libstdcxx_printers(gdb.current_objfile())

import re
import numpy

def key_str(key):
    rank_digits = 8
    id_digits = 56
    id_mask = 0x00ffffffffffffff

    types = ("NODE","LINE","SURFACE","ELEMENT","PARTICLE","CONSTRAINT")

    keyid = key & id_mask
    rank = key >> id_digits

    return '%s[%d]' % (types[rank],keyid)
    

class RelationPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        entity = self.val['m_entity'].dereference()
        key = int(entity['m_entityImpl']['m_key']['key'])
        rid = int(self.val['m_attr']['normal_view']['identifier'])
        rank = int(self.val['m_attr']['normal_view']['entity_rank'])
        return "(%s,%d,%d)" % (key_str(key),rid,rank)
        
    def display_hint(self):
        return 'Relation'

class EntityKeyPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        key = int(self.val['key'])
        return key_str(key)
        
    def display_hint(self):
        return 'EntityKey'
    
class BucketPtrPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        bucket = self.val.dereference()
        bucket = bucket['m_bucketImpl']
        if bucket['m_capacity']==0:
            return "bucket_nil"
        key = bucket['m_key']
        count = int(key[0])
        parts = [int(key[i+1]) for i in range(count)]
        #return "Bucket(%s)" % (",".join([str(i) for i in parts]),)
        mesh = bucket['m_mesh']
        meta = mesh['m_mesh_meta_data']
        universal = meta['m_universal_part'].dereference()
        allparts = universal['m_partImpl']['m_subsets']['_M_impl']['_M_start']
        myparts = [allparts[i].dereference() for i in parts[:-1]]
        partstr = ",".join([str(p['m_partImpl']['m_name']) for p in myparts])
        return "Bucket(%s,%d)" % (partstr,parts[-1])
    
    def display_hint(self):
        return 'Bucket'
    
class PartPtrPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        part = self.val.dereference()
        name = part['m_partImpl']['m_name']
        ordinal = part['m_partImpl']['m_universe_ordinal']
        return "Part(%s,%d)" % (name,ordinal)
    
    def display_hint(self):
        return 'Part'
    
class EntityPtrPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        if self.val==0:
            return "NULL"
        entity = self.val.dereference()
        key = int(entity['m_entityImpl']['m_key']['key'])
        return key_str(key)
    
    def display_hint(self):
        return 'Entity'


class EntityProcPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        proc = self.val['second']
        if self.val['first']==0:
            return '(NULL, %d)' % (proc)
        entity = self.val['first'].dereference()
        key = int(entity['m_entityImpl']['m_key']['key'])
        return '(%s, %d)' % (key_str(key),proc)
    
    def display_hint(self):
        return 'EntityProc'


class EntityKeyPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        try:
            # this handles EntityKey unions
            key = int(self.val['key'])
        except:
            # or we might have plain numbers
            key = int(self.val)
        return 'key ' + key_str(key)
    
    def display_hint(self):
        return 'EntityKey'


class KeyProcPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        key = int(self.val['first'])
        proc = self.val['second']
        return '(key %s, %d)' % (key_str(key),proc)
    
    def display_hint(self):
        return 'KeyProc'

class StdStringPrinter:
    "Print a std::basic_string of some kind"

    def __init__(self, val):
        self.val = val

    def to_string(self):
        # Make sure &string works, too.
        type = self.val.type
        if type.code == gdb.TYPE_CODE_REF:
            type = type.target ()

        # Calculate the length of the string so that to_string returns
        # the string according to length, not according to first null
        # encountered.
        ptr = self.val ['_M_dataplus']['_M_p']
        realtype = type.unqualified ().strip_typedefs ()
        reptype = gdb.lookup_type (str (realtype) + '::_Rep').pointer ()
        header = ptr.cast(reptype) - 1
        len = header.dereference ()['_M_length']
        #return self.val['_M_dataplus']['_M_p'].lazy_string (length = len)
        return self.val['_M_dataplus']['_M_p'].string ()

    def display_hint (self):
        return 'string'


class Epetra_MapPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        data = self.val['BlockMapData_'].dereference()
        ge = data['MyGlobalElements_']
        a = ge['A_']
        m = int(ge['M_'])
        gids = [str(a[i]) for i in xrange(m)]
        return " ".join(gids)

    def display_hint (self):
        return 'Epetra_Map'

class Epetra_IntSerialDenseVectorPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        a = self.val['A_']
        m = int(self.val['M_'])
        n = int(self.val['N_'])
        data = [[int(a[i+m*j]) for j in xrange(n)] for i in xrange(m)]
        array = numpy.array(data)
        return str(array)

    def display_hint (self):
        return 'Epetra_IntSerialDenseVector'

class Epetra_SerialDenseMatrixPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        a = self.val['A_']
        m = int(self.val['M_'])
        n = int(self.val['N_'])
        data = [[float(a[i+m*j]) for j in xrange(n)] for i in xrange(m)]
        array = numpy.array(data)
        return str(array)

    def display_hint (self):
        return 'Epetra_SerialDenseMatrix'


def point(p):
    x = p['x_']
    x = (float(x[0]),float(x[1]),float(x[2]),)
    return "(%.20g,%.20g,%.20g)" % x

def node(n):
    p = n['point_']
    return point(p)
    
class PointPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        p = self.val.dereference()
        cut_edges = p['cut_edges_']
        cut_sides = p['cut_sides_']
        return point(p)
        #return '%s\n\tedges={%s}\n\tsides={%s}' % (point(p),cut_edges,cut_sides)
    
    def display_hint(self):
        return 'Point'

class NodePrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        n = self.val.dereference()
        return node(n)
    
    def display_hint(self):
        return 'Node'

class EdgePrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        e = self.val.dereference()
        nodes = e['nodes_']
        return nodes
    
    def display_hint(self):
        return 'Edge'

class SidePrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        s = self.val.dereference()
        nodes = s['nodes_']
        return nodes
    
    def display_hint(self):
        return 'Side'

class LinePrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        s = self.val.dereference()
        p1 = s['p1_']
        p2 = s['p2_']
        cut_sides = s['cut_sides_']
        cut_elements = s['cut_elements_']
        return "L[%s-%s %s %s]" % (point(p1),point(p2),cut_sides,cut_elements)
    
    def display_hint(self):
        return 'Line'

class ElementPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        e = self.val.dereference()
        eid = e['eid_']
        cells = e['cells_']
        #nodes = e['nodes_']
        #nodes = "-".join([node(n) for n in nodes])
        #return "E%d[%s]" % (eid,nodes)
        #return "E%d " % eid
        return cells
    
    def display_hint(self):
        return 'Element'

class FacetPrinter:

    def __init__(self, val):
        self.val = val
        
    def to_string(self):
        s = self.val.dereference()
        planar = s['planar_']
        #if planar:
        points = s['points_']
        #else:
        #    points = s['triangulation_']
        return points
    
    def display_hint(self):
        return 'Facet'

class VolumeCellPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        v = self.val.dereference()
        facets = v['facets_']
        return facets
    
    def display_hint(self):
        return 'VolumeCell'

class LineSegmentPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        ls = self.val.dereference()
        points = ls['facet_points_']
        return points
    
    def display_hint(self):
        return 'LineSegment'

class TetHandlePrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        handle = self.val.dereference()
        points = handle['points_']
        return points
    
    def display_hint(self):
        return 'TetHandle'

class TetEntityPrinter:

    def __init__(self, val):
        self.val = val

    def to_string(self):
        entity = self.val.dereference()
        points = entity['handle_']['points_']
        return points
    
    def display_hint(self):
        return 'TetEntity'

def register_fe_printers (obj):
    "Register fe pretty-printers with objfile Obj."

    if obj == None:
        obj = gdb

    obj.pretty_printers.insert (0,lookup_function)


def lookup_function (val):
    "Look-up and return a pretty-printer that can print val."

    # Get the type.
    type = val.type

    # If it points to a reference, get the reference.
    if type.code == gdb.TYPE_CODE_REF:
        type = type.target ()

    # Get the unqualified type, stripped of typedefs.
    type = type.unqualified ().strip_typedefs ()

    # Get the type name.    
    typename = type.tag
    if typename == None:
        typename = str(type)

    #print typename

    # Iterate over local dictionary of types to determine
    # if a printer is registered for that type.  Return an
    # instantiation of the printer if found.
    for function in pretty_printers_dict:
        if function.search (typename):
            return pretty_printers_dict[function] (val)
      
    # Cannot find a pretty printer.  Return None.
    return None

pretty_printers_dict = {}

# any 'unsigned long long' is supposed to be an EntityKey
pretty_printers_dict[re.compile('^unsigned long long$')] = EntityKeyPrinter
pretty_printers_dict[re.compile('^std::pair<unsigned long long, int>$')] = KeyProcPrinter

pretty_printers_dict[re.compile('^stk::mesh::Entity \\*$')] = EntityPtrPrinter
pretty_printers_dict[re.compile('^const stk::mesh::Entity \\*$')] = EntityPtrPrinter
pretty_printers_dict[re.compile('^std::pair<stk::mesh::Entity\\*, unsigned int>$')] = EntityProcPrinter

pretty_printers_dict[re.compile('^stk::mesh::Bucket \\*$')] = BucketPtrPrinter
pretty_printers_dict[re.compile('^stk::mesh::Part \\*$')] = PartPtrPrinter

pretty_printers_dict[re.compile('^stk::mesh::Relation$')] = RelationPrinter
pretty_printers_dict[re.compile('^stk::mesh::EntityKey$')] = EntityKeyPrinter

# this is a modified copy
# needed since my gdb seems a little older than the pretty printers
pretty_printers_dict[re.compile('^std::basic_string<.*>$')] = lambda val: StdStringPrinter(val)

pretty_printers_dict[re.compile('^Epetra_Map$')] = Epetra_MapPrinter
pretty_printers_dict[re.compile('^Epetra_IntSerialDenseVector$')] = Epetra_IntSerialDenseVectorPrinter
pretty_printers_dict[re.compile('^Epetra_SerialDenseMatrix$')] = Epetra_SerialDenseMatrixPrinter

pretty_printers_dict[re.compile('^GEO::CUT::Point \\*$')] = PointPrinter
pretty_printers_dict[re.compile('^GEO::CUT::Node \\*$')]  = NodePrinter
pretty_printers_dict[re.compile('^GEO::CUT::Edge \\*$')]  = EdgePrinter
pretty_printers_dict[re.compile('^GEO::CUT::Side \\*$')]  = SidePrinter
pretty_printers_dict[re.compile('^GEO::CUT::Line \\*$')]  = LinePrinter
pretty_printers_dict[re.compile('^GEO::CUT::Element \\*$')]  = ElementPrinter
pretty_printers_dict[re.compile('^GEO::CUT::Facet \\*$')]  = FacetPrinter
pretty_printers_dict[re.compile('^GEO::CUT::VolumeCell \\*')] = VolumeCellPrinter

pretty_printers_dict[re.compile('^GEO::CUT::LineSegment \\*')] = LineSegmentPrinter

pretty_printers_dict[re.compile('^GEO::CUT::TetMesh::Handle<2> \\*')] = TetHandlePrinter
pretty_printers_dict[re.compile('^GEO::CUT::TetMesh::Handle<3> \\*')] = TetHandlePrinter
pretty_printers_dict[re.compile('^GEO::CUT::TetMesh::Handle<4> \\*')] = TetHandlePrinter

pretty_printers_dict[re.compile('^GEO::CUT::TetMesh::Entity<2> \\*')] = TetEntityPrinter
pretty_printers_dict[re.compile('^GEO::CUT::TetMesh::Entity<3> \\*')] = TetEntityPrinter
pretty_printers_dict[re.compile('^GEO::CUT::TetMesh::Entity<4> \\*')] = TetEntityPrinter

register_fe_printers(gdb.current_objfile())
