import argparse
import ConfigParser
import os

from Bio.PDB import *
from PDBUtils import PDBUtils

from Sphere import Sphere


DEBUG=0 # debugging flag

if (__name__ == '__main__'):
  #
  # read command-line arguments
  argsparser = argparse.ArgumentParser(description='SphereGrinder command-line arguments') 
  argsparser.add_argument('--model', '-m', \
    help='path to file with model PDB structure', required=True)
  argsparser.add_argument('--reference', '-r',  \
    help="path to file with reference PDB structure", required=True)
  argsparser.add_argument('--output', '-o',  \
    help="path to output file", required=True)
  args = argsparser.parse_args()
  model = args.model
  target = args.reference
  outf = args.output
  #
  # read default config file
  config = ConfigParser.ConfigParser()
  try:
    config.read("default.ini")
    radii_str = config.get('Sphere_params', 'radii')
    caOnly = config.getboolean('Sphere_params', 'caOnly')
    incompleteBBmax = config.getin('PDB_structure', 'incompleteBBmax')
  except (OSError("default.ini"), KeyError, IOError, ConfigParser.NoSectionError): 
    # if config file doesn't exists, can't be read, 
    # doesn't have proper parameters, etc
    radii_str = '1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,300.0'
    caOnly = False # by default process all atoms
    incompleteBBmax = 5
  # convert string of parameter radii into array of floats
  radii = []
  for r in radii_str.split(","):
    try:
      radii.append(float(r))
      if (r <= 0.0):
        raise ValueError
    except ValueError:
        raise ValueError("Radius value {} is in wrong format. Positive float is expected".format(r))
  radii = sorted(radii)
  if DEBUG == 1:
    for r in radii:
      print "{:.1f},".format(r),
    print ""
    sys.exit(0)
  #
  # read structures
  mf = open(model, "r")
  parser = PDBParser()
  model = parser.get_structure('model', mf)
  mf.close()
  tf = open(target, "r")
  target = parser.get_structure('target', tf)
  tf.close()
  #
  #print model in DEBUG mode
  if DEBUG == 1:
    m_atoms = model.get_atoms()
    for a in m_atoms:
      print "Res:{} Atom:{} x:{} y:{} z:{}".format(a.get_parent().get_resname(), a.get_name(), a.get_coord()[0], a.get_coord()[1], a.get_coord()[2])
  # clean structures: 
  # - remove non-standard atoms
  # - extract only first chain from first model 
  # - check completeness of backbones
  # - check alignment of model and reference 
  pdbu = PDBUtils()
  model = pdbu.filterNonStdAtoms(model)
  target = pdbu.filterNonStdAtoms(target)
  model = pdbu.checkBB(model, incompleteBBmax)
  target = pdbu.checkBB(target)
  m_ch = pdbu.extractFirstChain(model)
  t_ch = pdbu.extractFirstChain(target)
  if DEBUG==1:
    m_atoms = model.get_atoms()
    print "\n\n\nAFTER FILTERING!!!\n\n\n"
    for a in m_atoms:
      print "Res:{} Atom:{} x:{} y:{} z:{}".format(a.get_parent().get_resname(), a.get_name(), a.get_coord()[0], a.get_coord()[1], a.get_coord()[2])
    t_residues = target.get_residues() 
    i = 0
    for t_res in t_residues:
      rand_res = t_res
      if i == 10:
       break
      i = i + 1
    print "!!! !!!"
    sph = Sphere(m_ch, t_ch)
    pr = sph.processResRadius(rand_res, 300.0)
    print "Before swap: rmsd:{}|frac:{}".format(pr[0], pr[1])
    sph.swapAmbiguousAtomLabels()
    pr = sph.processResRadius(rand_res, 300.0)
    print "After swap: rmsd:{}|frac:{}".format(pr[0], pr[1])
  #
  # SphereGrinder engine
  sph = Sphere(m_ch, t_ch)
  sph.swapAmbiguousAtomLabels()
  data = sph.process(radii, caOnly)
  #
  # output results to file
  outdir = os.path.dirname(outf)
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  if not os.path.isdir(outdir):
    raise OSError("Can't make dir {}".format(outdir))
  try:
    outh = open(outf, "w+")
    outh.write(sph.toString(data))
    outh.write("\n\nGlobal score: {:.5f}".format(sph.calcGlobalScore(data)))
    outh.close()
  except OSError:
    raise OSError("Failed to write to file {}".format(outf))
  
