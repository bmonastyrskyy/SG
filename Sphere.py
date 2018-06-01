from Bio.PDB import *
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
from Bio.PDB import Superimposer
from PDBUtils import PDBUtils

class Sphere:

  def __init__(self, model_ch, ref_ch):
    """ The initialiser of the class.
	
	Arguments:
	 arg1 - model chain, instance of class Bio.PDB.Chain
	 arg2 - reference chain, instance of class Bio.PDB.Chain
        
     The method initiates an instance of the class with two attributes:
      - model_chain
      - ref_chain 
    """
    if (model_ch is None ):
      raise TypeError("The initialiser of the class Sphere expects two non-None arguments of class Bio.PDB.Chain")
    if (not model_ch.__class__.__name__ == 'Chain'):
      raise TypeError("The initialiser of the class Sphere expects the first argument of class Bio.PDB.Chain")
    self.model_chain = model_ch
    if (ref_ch is None ):
      raise TypeError("The initialiser of the class Sphere expects two non-None arguments of class Bio.PDB.Chain")
    if (not ref_ch.__class__.__name__ == 'Chain'):
      raise TypeError("The initialiser of the class Sphere expects the second argument of class Bio.PDB.Chain")
    self.ref_chain = ref_ch


  
  def getReferenceAtomsWithinSphere(self, residue, radius, caOnly=False):
    if (not residue.__class__.__name__ == 'Residue'):
      raise TypeError("The function process(residue, radius, caOnly) expects the first argument of class Bio.PDB.Residue")
    ns = NeighborSearch(Selection.unfold_entities(self.ref_chain, 'A'))
    ref_atoms = ns.search(residue['CA'].get_coord(), radius, level='A')
    if (caOnly):
      cas = []
      for a in ref_atoms:
        if (a.get_name() == 'CA'):
           cas.append(a)
      return cas
    return ref_atoms


  
  def getCorrespondAtomsInModel(self, ref_atoms):
    common_atoms_ref = []
    common_atoms_model = []
    for ref_a in ref_atoms:
      res = ref_a.get_parent()
      res_id = res.get_id()[1] # resseq
      atom_id = ref_a.get_id() # atom name
      try:
        model_a = self.model_chain[res_id][atom_id]
        if ('swap' in model_a.xtra):
          if (model_a.xtra['swap'] == True):
            atom_id = PDBUtils.ALTERNATIVE_ATOMS[res.get_resname()][atom_id]
            model_a = self.model_chain[res_id][atom_id]
        common_atoms_ref.append(ref_a)
        common_atoms_model.append(model_a)
      except KeyError: # if atom doesn't exist in model - skip it
        continue
    return (common_atoms_ref, common_atoms_model, (len(common_atoms_ref)*1.0)/(1.0*len(ref_atoms)))


  
  def calcRMSD(self, fixed_atoms, moving_atoms):
    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    return sup.rms


  
  def processResRadius(self, residue, radius, caOnly=False):
    """ Method 
	- identifies atoms within a sphere with center at CA of the residue and given radius
	- picks up the corresponding atoms in model
	- calculates RMSD and fraction of atoms
	
	Arguments:
	residue - where the sphere centered, of class Bio.PDB.Residue
	radius - self-explonatory, of type float
        caOnly - boolean flag whether to pick up all atoms or Ca's only
    """
    if (not residue.__class__.__name__ == 'Residue'):
      raise TypeError("The function processresRadius(residue, radius, caOnly) expects the first argument of class Bio.PDB.Residue")
    try: # if missing residue in model
      self.model_chain[residue.get_id()[1]]
    except KeyError:
      return (None, None)
    ref_atoms = self.getReferenceAtomsWithinSphere(residue, radius, caOnly)
    corresp = self.getCorrespondAtomsInModel(ref_atoms)
    if (len(corresp[0]) == 0):
      return (None, None);
    rms = self.calcRMSD(corresp[0], corresp[1])
    return (rms, corresp[2])



  def processRes(self, residue, radii, caOnly=False):
    """ 
       radii - list of radii to process
       The method returns dictionary of tulips: key - radius, value - (rms, fraction)
    """
    if (not residue.__class__.__name__ == 'Residue'):
      raise TypeError("The function processRes(residue, radii, caOnly) expects the first argument of class Bio.PDB.Residue")
    result = [None]*len(radii)
    for index,radius in enumerate(radii):
      radius = float("{:.1f}".format(radius))
      (rms,frac) = self.processResRadius(residue, radius, caOnly)
      result[index] = {'radius':radius, 'rms': rms, 'frac':frac}
    return result
  
  
  
  def process(self, radii, caOnly=False):
    """
      The method returns the dictionary of the following structure: 
      {
         resseq : {
           'resname': 3LETTER_AA_CODE,
           'sph': [
               {'radius': radius,
                    'rms': rms,
                    'frac': frac
               }
          ]
         } 
      }
    """
    result = {} 
    for residue in self.ref_chain: # loop over residues
      resseq = residue.id[1]  # residue sequence number 
      resname = residue.get_resname() # residue name
      sphs = self.processRes(residue, radii, caOnly) 
      result[resseq] = {'resname':resname,'sph':sphs}
    return result


  def toString(self, data):
    """
      The method formats string from dictionary - output of process.
      Expected argument - dictionary - output of the method process(radii, caOnly)
    """
    # header flag
    header = True
    result = ''
    for resseq in sorted(data.iterkeys()):
      if (header):
        result = "NAME"
        for index in xrange(len(data[resseq]['sph'])):
          result = result + ";" + "sphere-{:.1f}".format(data[resseq]['sph'][index]['radius'])
        result = result + "\n"
        header= False
      result = result + "{}\t{}".format(resseq,data[resseq]['resname'])
      for index in xrange(len(data[resseq]['sph'])):
        try:
          result = result + ";" + "{:.8f}|{:.4f}".format(\
          data[resseq]['sph'][index]['rms'], \
          data[resseq]['sph'][index]['frac']\
          )
        except ValueError:
          result = result + ";" + "NN|NN"
      result = result + "\n"
    return result
  
  def calcGlobalScore(self, data, radius_index=5, thresholds=[2.0, 4.0]):
    """
    The method calculates the global score:
    average_over_thresholds(sum of fractions above the rmsd_threshold normilized by the length of the reference structure)
    
    : param data
    : type dictionary : output of method process()

    : param radius_index : the index of the element of the list "radii"
    : type int : default value 5 corresponds to the sphere radius 6.0 
                 for default value of radii = '1.0,2.0,...,30.0,300.0'

    : param thresholds : rmsd threshods
    : type list : default [2.0, 4.0]
    """
    # loop over residues in reference structure
    n=0 # counter of residues in reference structure
    tmp = [0.0 for i in xrange(len(thresholds))] # auxilary list to store temporary results
    for resseq in sorted(data.iterkeys()):
      n = n + 1
      for i,thresh in enumerate(thresholds):
        try:
          if data[resseq]['sph'][radius_index]['rms']<=thresh:
            tmp[i] = tmp[i] + data[resseq]['sph'][radius_index]['frac']
        except ValueError: 
          # do nothing if the vaue not a number 
          # that might be the case of missing residues/atoms in model
          pass
    # normalized by length and calculate average
    try:
      return float(sum(tmp))/(n*len(tmp))
    except ZeroDivisionError:
      return 0.0
  
  def swapAmbiguousAtomLabels(self):
    """
    The method checks if the atom with alternative label should be swapped/relabeled.
    The decision is made based on comparison of the rmsd for two sets of atoms:
        - all Ca's + original labels for the pair of atoms in question;
        - all Ca's + swapped labels for the pair of atoms;
    The lower rmsd makes a winner.
    """
    cas_ref = []
    cas_model = []
    atoms2check = []
    # loop over all atoms in reference structure
    # select CA's and identify atoms to be checked (with ambiguious labels)
    for ref_a in Selection.unfold_entities(self.ref_chain, 'A'):
      res = ref_a.get_parent()
      res_id = res.get_id()[1] # resseq
      atom_id = ref_a.get_id() # atom name
      if atom_id == 'CA':  # select CA only
        try:
          model_a = self.model_chain[res_id][atom_id]
          cas_ref.append(ref_a)
          cas_model.append(model_a)
        except KeyError: # if atom doesn't exist in model - skip it
          continue
      elif res.get_resname() in PDBUtils.ALTERNATIVE_ATOMS: # if residue 'ASP' || 'GLU' || 'TYR' || 'PHE' || 'ARG'
        try:
          a1 = atom_id # atom_name
          a2 = PDBUtils.ALTERNATIVE_ATOMS[res.get_resname()][a1] # alternative atom_name, e.g. OD2 for OD1 in ASP
          # check if all aleternative atom are present in reference and model 
          atoms2check.append(ref_a)
          #print "{}".format(ref_a.get_name())
        except KeyError:
          continue
    # loop over all atoms to be checked
    for a2c in atoms2check:
      res = a2c.get_parent()
      res_id = res.get_id()[1] # resseq
      atom_id = a2c.get_id()
      a1 = a2c.get_id()
      a2 = PDBUtils.ALTERNATIVE_ATOMS[res.get_resname()][a1]
      try:
        self.model_chain[res_id]
      except KeyError:
        continue
      if (a1 in self.model_chain[res_id]):
        if ('swap' in self.model_chain[res_id][a1].xtra): # we have processed this case
          continue
      if (a2 in self.model_chain[res_id]):
        if ('swap' in self.model_chain[res_id][a2].xtra): # we have processed this case
          continue
      if ((a2 in self.ref_chain[res_id]) and (a1 in self.model_chain[res_id]) and (a2 in self.model_chain[res_id])):
        atoms_ref = list(cas_ref)
        atoms_ref.append(self.ref_chain[res_id][a1])
        atoms_ref.append(self.ref_chain[res_id][a2])
        # orig labeling
        atoms_m1 = list(cas_model)
        atoms_m1.append(self.model_chain[res_id][a1])
        atoms_m1.append(self.model_chain[res_id][a2])
        rms1 = self.calcRMSD(atoms_ref, atoms_m1)
        # swapped labeling
        atoms_m2 = list(cas_model)
        atoms_m2.append(self.model_chain[res_id][a2])
        atoms_m2.append(self.model_chain[res_id][a1])
        rms2 = self.calcRMSD(atoms_ref, atoms_m2)
        if rms1 > rms2: # the labels should be swapped
          self.model_chain[res_id][a1].xtra['swap'] = True
          self.model_chain[res_id][a2].xtra['swap'] = True
        else:
          self.model_chain[res_id][a1].xtra['swap'] = False
          self.model_chain[res_id][a2].xtra['swap'] = False
      elif ((a2 in self.ref_chain[res_id]) and (a1 in self.model_chain[res_id])):
        atoms_ref1 = list(cas_ref)
        atoms_ref1.append(self.ref_chain[res_id][a1])
        atoms_ref2 = list(cas_ref)
        atoms_ref2.append(self.ref_chain[res_id][a2])
        # orig labeling
        atoms_m = list(cas_model)
        atoms_m.append(self.model_chain[res_id][a1])
        rms1 = self.calcRMSD(atoms_ref1, atoms_m)
        rms2 = self.calcRMSD(atoms_ref2, atoms_m)
        if rms1 > rms2: # the labels should be swapped
          self.model_chain[res_id][a1].xtra['swap'] = True
        else:
          self.model_chain[res_id][a1].xtra['swap'] = False
      elif ((a2 in self.ref_chain[res_id]) and (a2 in self.model_chain[res_id])):
        atoms_ref1 = list(cas_ref)
        atoms_ref1.append(self.ref_chain[res_id][a1])
        atoms_ref2 = list(cas_ref)
        atoms_ref2.append(self.ref_chain[res_id][a2])
        # orig labeling
        atoms_m = list(cas_model)
        atoms_m.append(self.model_chain[res_id][a2])
        rms1 = self.calcRMSD(atoms_ref1, atoms_m)
        rms2 = self.calcRMSD(atoms_ref2, atoms_m)
        if rms1 > rms2: # the labels should be swapped
          self.model_chain[res_id][a2].xtra['swap'] = True
        else:
          self.model_chain[res_id][a2].xtra['swap'] = False
      elif ((a1 in self.model_chain[res_id]) and (a2 in self.model_chain[res_id])):
        atoms_ref = list(cas_ref)
        atoms_ref.append(self.ref_chain[res_id][a1])
        # orig labeling
        atoms_m1 = list(cas_model)
        atoms_m1.append(self.model_chain[res_id][a1])
        rms1 = self.calcRMSD(atoms_ref, atoms_m1)
        # swapped labeling
        atoms_m2 = list(cas_model)
        atoms_m2.append(self.model_chain[res_id][a2])
        rms2 = self.calcRMSD(atoms_ref, atoms_m2)
        if rms1 > rms2: # the labels should be swapped
          self.model_chain[res_id][a1].xtra['swap'] = True
          self.model_chain[res_id][a2].xtra['swap'] = True
        else:
          self.model_chain[res_id][a1].xtra['swap'] = False
          self.model_chain[res_id][a2].xtra['swap'] = False
     
