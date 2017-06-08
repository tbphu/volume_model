import numpy as np
import libsbml

class OptimizationBounds(object):
    def __init__(self, xmax, xmin):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


def get_initial_assignments_dict(model):
	sbml = model.getSBML()
	sbml_document = libsbml.readSBMLFromString(sbml)
	sbml_model = sbml_document.getModel()
	assignments = {}
	for ia in sbml_model.getListOfInitialAssignments():
		assignments[ia.getId()] = libsbml.formulaToString(ia.getMath())
	return assignments