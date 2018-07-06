# user-defined Exceptions
class Error(Exception):
  """ Base class for other exceptions"""
  pass

class WrongAminoAcidError(Error):
  """ Raised when the amino acid at the certain position in the model isn't the same as that in the target """
  pass

class IncompleteBBError(Error):
  """ Raised when a amino acid has incomplete set of backbome atoms"""
  pass
