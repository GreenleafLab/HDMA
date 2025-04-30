# Author: Selin Jessa
# 2024 

import torch
from bpnetlite.attribute import deep_lift_shap

class ChromBPNetWrapper(torch.nn.Module):
	"""A wrapper class that returns counts and transformed profiles
	from a bias-corrected ChromBPNet model.

	This class takes in a trained model, specifically expecting a BPNet model
	from a trained bias-corrected ChromBPNet model (i.e. `chrombpnet_nobias.h5`),
	_not_ the full accessibility model.

	It returns outputs the same shape as the BPNet model, but with the predicted
	profile logits softmaxed and scaled by the exponent of the predicted counts.
	This is for convenience when plotting, since in e.g. bigwigs, we are inspecting
	the scaled profile, not the unnormalized logits.

	Parameters
	----------
	model: torch.nn.Module
		A torch model to be wrapped.
	"""

	def __init__(self, model):
		super().__init__()
		self.model = model
	
	def forward(self, X):
        
		y = self.model(X)

		# predicted profile logits
		y_profile = y[0]

		# predicted log counts
		y_counts = y[1]
              
		# softmax the logits to get profile probabilities
		y_profile = torch.nn.functional.softmax(y_profile, dim = -1)

		# scale the profile by the exponentiated predicted ln(counts)
		y_profile = y_profile * torch.exp(y_counts).unsqueeze(2)

		return y_profile, y_counts
	



def reshape_folds(y):
	"""
	Separates a list of [profile, counts] lists, one for each of the model folds, into two separate lists:
	one for profiles and one for counts, with each element corresponding to a fold, then stacks the outputs.

	Args:
		y (list): A list where each element is a list containing two items:
		[profile (torch.Tensor or similar), counts (torch.Tensor or similar)].

	Returns:
		tuple: A tuple containing two lists:
				(profiles (list), counts (list))
				profiles: A list of the profile elements from the input data.
				counts: A list of the counts elements from the input data.
	"""
	profiles = []
	counts = []
	
	for item in y:
		if len(item) == 2:  # Ensure each sublist has two elements
			profiles.append(item[0])  # Assuming profile is the first element
			counts.append(item[1])  # Assuming counts is the second element
		else:
			print("Warning: Sublist has an unexpected number of elements.")

	# stack the profiles and counts into a new dimension
	profiles = torch.stack(profiles, dim=0)
	counts = torch.stack(counts, dim=0)
	
	return profiles, counts
