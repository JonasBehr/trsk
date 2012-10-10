"""
simple script to train splice site predictors
"""

from arts2_factory import create_hashed_features_wdk
from shogun_factory_new import create_labels, create_svm
from helper import Options, assess
from io_pickle import save, load



def parse_preprocessing_output(input_fn):
	"""
	parses fasta-oid format:

	%-1 aaact
	%1 ttact
	%-1 aaact

	@param input_fn: path to input file
	"""
	input_file = file(input_fn)

	examples = []
	labels = []
	first = 1

	for line in input_file:
		if line.startswith("#"):
			tokens = line[1:].strip().split(" ")
			if len(tokens) == 2:	
				if first == 1:
					exlen = len(tokens[1]); 
					first = 0
					print exlen
				alphabet = set(tokens[1])
				should_be_alphabet = set(["a", "g", "t", "c", "n", "N"])
				if not alphabet.issubset(should_be_alphabet):
		   			print "ARGH", alphabet.symmetric_difference(should_be_alphabet)
				elif not len(tokens[1]) == exlen: 
					print line
				else: 
					examples.append(tokens[1].upper().replace("N", "A"))
		   			labels.append(int(tokens[0]))
			else: 
				print len(tokens), line

	return (examples, labels)	



def train_splice_predictor(examples, labels, param):
	"""

	@param examples: list of strings
	@param labels: list of integers {-1,1}
	"""
 
	
	##########################
	#   build classifier
	##########################

	feat_train = create_hashed_features_wdk(param.flags, examples)
	lab = create_labels(labels)
	svm = create_svm(param, feat_train, lab)

	svm.train()

	return svm

def init_predictor(examples, labels, param, w):
	"""

	@param examples: list of strings
	@param labels: list of integers {-1,1}
	@param w: weight vector of trained svm
	"""
 
	
	##########################
	#   build classifier
	##########################

	feat_train = create_hashed_features_wdk(param.flags, examples)
	lab = create_labels(labels)
	svm = create_svm(param, feat_train, lab)

	svm.set_w(w)

	return svm


def predict(predictor, examples, param):
	"""
	make prediction on examples using trained predictor

	@param predictor: trained predictor
	@type predictor: SVM object
	@param examples: list of examples
	@type examples: list 
	"""
	
	#shogun data
	feat = create_hashed_features_wdk(param.flags, examples)
	
	#predict
	svm_out = predictor.classify(feat).get_labels()

	return svm_out

def define_param():
	# flags
	flags = {}
	flags["normalize_cost"] = False
	flags["kernel_cache"] = 200
	flags["use_bias"] = False
	flags["epsilon"] = 0.01

	# arts params
	flags["svm_type"] = "liblineardual"
	flags["degree"] = 24
 

	#create mock param object by freezable struct
	param = Options()
	param.cost = 1.0
	param.id = 666
	param.flags = flags
	
	param.freeze()
	
	return param

def Main():
	"""
	main procedure
	"""

	##########################
	#   set parameters
	##########################

	param = define_param();

	(examples, labels) = parse_preprocessing_output("/fml/ag-raetsch/home/jonas/svn/tools/ngs/train.txt")
	
	trained_predictor = train_splice_predictor(examples, labels, param)

	w = trained_predictor.get_w()

	save("classifier_w.pickle", w)	

	w2 = load("classifier_w.pickle")

	predictor2 = init_predictor(examples, labels, param, w2)	

	(examples_test, labels_test) = parse_preprocessing_output("/fml/ag-raetsch/home/jonas/svn/tools/ngs/test.txt")


	out = predict(trained_predictor, examples_test, param)
	out2 = predict(predictor2, examples_test, param)

	target = "auPRC"
	performance = assess(out, labels_test, target)
	print "auPRC", performance

	target = "auROC"
	performance = assess(out, labels_test, target)
	print "auROC", performance

	target = "auPRC"
	performance = assess(out2, labels_test, target)
	print "auPRC", performance

	target = "auROC"
	performance = assess(out2, labels_test, target)
	print "auROC", performance


if __name__ == "__main__":
	Main()
