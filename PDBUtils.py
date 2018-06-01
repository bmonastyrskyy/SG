from Bio.PDB import *
from SGExceptions import *

class PDBUtils:
  """
  Dictionary of standard heavy atoms: key - 3letter aa_code, value - tuple of atom names
  """
  AA_STANDARD_ATOMS = {
    'ALA': ('N', 'CA', 'C', 'O', 'CB'),
    'ARG': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'),
    'ASN': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'),
    'ASP': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'),
    'CYS': ('N', 'CA', 'C', 'O', 'CB', 'SG'),
    'GLU': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'),
    'GLN': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'),
    'GLY': ('N', 'CA', 'C', 'O'),
    'HIS': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'),
    'ILE': ('N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'),
    'LEU': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'),
    'LYS': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'),
    'MET': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'),
    'PHE': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
    'PRO': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'),
    'SER': ('N', 'CA', 'C', 'O', 'CB', 'OG'),
    'THR': ('N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'),
    'TRP': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'),
    'TYR': ('N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'),
    'VAL': ('N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2')
  }
  """
  Dictionary of alternative atoms - those which labels can be switched.
  """
  ALTERNATIVE_ATOMS = {
    'ARG': {
             'NH1': 'NH2',
             'NH2': 'NH1'
           },
    'ASP': {
             'OD1': 'OD2',
             'OD2': 'OD1'
           },
    'GLU': {
             'OE1': 'OE2',
             'OE2': 'OE1'
           },
    'LEU': {
             'CD1': 'CD2',
             'CD2': 'CD1'
           },
    'TYR': {
             'CD1': 'CD2',
             'CD2': 'CD1',
             'CE1': 'CE2',
             'CE2': 'CE1'
           },
    'PHE': {
             'CD1': 'CD2',
             'CD2': 'CD1',
             'CE1': 'CE2',
             'CE2': 'CE1'
           },
    'VAL': {
             'CG1': 'CG2',
             'CG2': 'CG1'
           }
  }
  """
  Dictionary maps 3LETTER code to 1LETTER code
  """
  A321 = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
  }
  """
  Dictionary maps 1LETTER code to 3LETTER code
  """
  A123 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'E': 'GLU',
    'Q': 'GLN',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL'
  }
  """
  Tuple of backbome atoms
  """
  BB = ('N', 'C', 'CA', 'O')
  
  def __init__(self):
    """
    The initialiser of the class.
    There are no attributes.
    """
    pass
  
  def filterNonStdAtoms(self, struct):
    """ The method filters out all non-standard atoms from all residues.
       
    """
    for m in struct: # loop over pdb models
      for ch in m:   # loop over chains
        for res in list(ch):  # loop over residues
          # if residue is not standard - remove it
          if not res.get_resname() in self.AA_STANDARD_ATOMS:
            ch.detach_child(res.id)
            continue
          for a in list(res):  # loop over atoms
            try:
              if not a.get_id() in self.AA_STANDARD_ATOMS[res.get_resname()] :
                res.detach_child(a.get_id())
                #print "Res:{} Atom:{}".format(res.get_resname(), a.get_id())
            except Exception:
                pass
    return struct;
  
  def extractFirstChain(self, struct):
    """ The method extacts first chain from the first model from Bio.PDB.Structure """
    for ch in struct[0]:
      return ch

  def checkAlignTwoChains(self, m_ch, ref_ch):
    """ The  method checks the alignment of two chains: prediction and referenece.
	The residues have to be matched at the corresponding posistions,
	otherwise the WrongAminoAcidException is raised and the program exits."""
    for res in ref_ch:
      if (not m_ch[res.get_id()[1]].get_resname() == res.get_resname()): 
        message = "ERROR! Wrong amino acid: in model {} {}, in target {} {}".format(\
        res.get_id()[1], m_ch[res.get_id()[1]].get_resname(), res.get_id()[1], \
        res.get_resname())
        raise WrongAminoAcidError(message)
  
  def checkBB(self, struct):
    """ Method checks if all residues have all backbone atoms. 
        Otherwise it raises IncompleteBBError.
    """
    for m in struct: # loop over pdb models
      for ch in m:   # loop over chains
        for res in ch:  # loop over residues
          for bb in self.BB:  # loop over atoms of backbone N, C, CA, O
            try:
              res[bb]
            except KeyError:
              message = "In residue {} {} backbone atom {} is missing.".format(\
              res.get_resname(), res.id[1], bb)
              raise IncompleteBBError(message)
    return struct;
