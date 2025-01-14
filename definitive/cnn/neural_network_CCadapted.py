# -*- coding: utf-8 -*-

"""
This code has been adapted from the original script by Nielsen in 2015.
Moreover, it has been further modified and adjusted to the coiled-coils data.
REFERENCES:
- http://neuralnetworksanddeeplearning.com/chap1.html
- https://github.com/mnielsen/neural-networks-and-deep-learning/blob/master/data/mnist.pkl.gz
- https://github.com/nchen909/NN-by-hand
"""

## CODE TO CONFIGURE THE NEURAL NETWORK:
# LIBRARIES IMPORTATION
from google.colab import drive
drive.mount('/content/drive')
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split

# FUNCTIONS

# To calculate the sigmoid function of the NN
def sigmoid(z):
    return 1.0/(1.0+np.exp(-z))

# To calculate the derivative of the sigmoid function, used for backpropagation purposes
def sigmoid_prime(z):
    return sigmoid(z)*(1-sigmoid(z))

# NEURAL NETWORK CREATION
# The NN is defined as an class-type object
class NeuralNetwork(object):
  def __init__(self, sizes): # Function to define the NN architecture
    self.layers = len(sizes) # Number of layers
    self.sizes = sizes # Number of nodes/neurons per layer
    # Omitting the first layer (input layer), the biases/weigths for the other hidden layers
    # are calculated considering a Gaussian distribution (mean:0, vaariance:1)
    # y will refer to the neurons number per layer; x to the number of weights per neuron (which is equal to
    # the number of inputs that neuron will use to compute the output)
    self.biases = [np.random.randn(y, 1) for y in sizes[1:]]
    self.weights = [np.random.randn(y, x) for x, y in zip(sizes[:-1], sizes[1:])]

  # The NN is a feedforward type network, defined by this function
  def feedforward(self, a):
    for bias, weight in zip(self.biases, self.weights): # Iteration though each layer's biases and weigths
      a = sigmoid(np.dot(weight, a)+bias)
    return a

  # The NN has a mini-batch stochastic gradient descent
  def SGD(self, training_data, epochs, mini_batch_size, eta,
          test_data=None):
      accuracy_list = []
      loss_list = []
      if test_data:
        n_test = len(test_data)
        n = len(training_data)
        for j in range(epochs):
            loss = 0
            epoch_loss = 0
            random.shuffle(training_data)
            mini_batches = [training_data[k:k+mini_batch_size]
                for k in range(0, n, mini_batch_size)]
            for mini_batch in mini_batches:
                self.update_mini_batch(mini_batch, eta)
            # Calculate loss for the whole training set (average loss over the epoch)
            for x, y in training_data:
                # Feedforward to get the output
                output = self.feedforward(x)
                # Compute the loss for this sample (using cost_derivative, which calculates the error)
                loss = np.sum(np.square(output - y))  # Sum of squared errors (this is one way to compute loss)
                epoch_loss += loss
            # Average loss for this epoch
            avg_loss = epoch_loss / n
            loss_list.append(avg_loss)  # Append the loss for this epoch only once per epoch

            if test_data:
                print(f"Epoch {j+1}: {self.evaluate(test_data)} / {n_test}. Accuracy: {round((self.evaluate(test_data))/n_test,2)}.")
                accuracy_list.append(round((self.evaluate(test_data))/n_test,2))
            else:
                print("Epoch {0} complete".format(j+1))
      return epochs, accuracy_list, loss_list

  # Updating the NN's biases and weights after applying SGD with backpropagation to one mini batch
  def update_mini_batch(self, mini_batch, eta):
        nabla_b = [np.zeros(b.shape) for b in self.biases]
        nabla_w = [np.zeros(w.shape) for w in self.weights]
        for x, y in mini_batch:
            delta_nabla_b, delta_nabla_w = self.backprop(x, y)
            nabla_b = [nb+dnb for nb, dnb in zip(nabla_b, delta_nabla_b)]
            nabla_w = [nw+dnw for nw, dnw in zip(nabla_w, delta_nabla_w)]
        self.weights = [w-(eta/len(mini_batch))*nw
                        for w, nw in zip(self.weights, nabla_w)]
        self.biases = [b-(eta/len(mini_batch))*nb
                       for b, nb in zip(self.biases, nabla_b)]

  # Application of backpropagation in the NN
  def backprop(self, x, y):
        nabla_b = [np.zeros(b.shape) for b in self.biases]
        nabla_w = [np.zeros(w.shape) for w in self.weights]
        # feedforward
        activation = x
        activations = [x] # list to store all the activations, layer by layer
        zs = [] # list to store all the z vectors, layer by layer
        for b, w in zip(self.biases, self.weights):
            z = np.dot(w, activation)+b
            zs.append(z)
            activation = sigmoid(z)
            activations.append(activation)
        # backward pass
        delta = self.cost_derivative(activations[-1], y) * \
            sigmoid_prime(zs[-1])
        nabla_b[-1] = delta
        nabla_w[-1] = np.dot(delta, activations[-2].transpose())
        for l in range(2, self.layers):
            z = zs[-l]
            sp = sigmoid_prime(z)
            delta = np.dot(self.weights[-l+1].transpose(), delta) * sp
            nabla_b[-l] = delta
            nabla_w[-l] = np.dot(delta, activations[-l-1].transpose())
        return (nabla_b, nabla_w)

  # Evaluate the correct predictions when test data is analyzed: test inputs
  # for which the NN has correctly predicted the result are returned.
  def evaluate(self, test_data):
        test_results = [(np.argmax(self.feedforward(x)), y)
                        for (x, y) in test_data]
        return sum(int(x == y) for (x, y) in test_results)

  # Finally, the vector of partial derivaties for the output activations is returned
  def cost_derivative(self, output_activations, y):
        return (output_activations-y)

# UPLOADING THE COILED-COILS DATASET
# Libraries
import numpy as np

def load_data():
    f = pd.read_csv('/content/drive/MyDrive/Semestre 3 (2024-2025)/lab_shared/data_overlap_prolinecorrected_toNN.csv',
                    sep="\t", header = None)
    data = f.to_numpy()
    # Specifically select number of digit 2
    data2 = data[(data[:,1]) == 2]
    data3 = data[(data[:,1]) == 3]
    data2_ind = np.random.choice(len(data2), size=188, replace=False)
    data2 = data2[data2_ind]
    # Specifically select digits 2 and 3
    #data = data[data2 | (data[:, 1] == 3)]
    # Join selected data:
    data = np.append(data2, data3, axis=0)
    print(len(data))
    training_data, temp_data = train_test_split(data, test_size=0.3, random_state=312)
    validation_data, test_data = train_test_split(temp_data, test_size=0.5, random_state=312)
    training_data = (training_data[:,0], training_data[:,1]-2) # Extract data and labels
    validation_data = (validation_data[:,0], validation_data[:,1]-2)
    test_data = (test_data[:,0], test_data[:,1]-2)
    """
    proportions_train = {
        2: 1,
        3: 1,
        4: 1,
        5: 1,
        6: 1,
        7: 1,
        8: 1
    }

    proportions_test = {
        2: 1,
        3: 1,
        4: 1,
        5: 1,
        6: 1,
        7: 1,
        8: 1
    }
    proportions_val = {
        2: 1,
        3: 1,
        4: 1,
        5: 1,
        6: 1,
        7: 1,
        8: 1
    }

    def sample_indices(data, proportions):
      # Get the class labels and data
      x, y = data
      # Initialize an empty list to store the indices for the selected data
      selected_indices = []

      # For each class, calculate the number of samples to select
      for class_label, proportion in proportions.items():
          class_indices = np.where(y == class_label)[0]
          num_samples = int(len(class_indices) * proportion)
          selected_class_indices = np.random.choice(class_indices, num_samples, replace=False)
          selected_indices.extend(selected_class_indices)
      # Return the selected indices
      return np.array(selected_indices)

    # Sample indices for training, validation, and test data
    train_indices = sample_indices(training_data, proportions_train)
    val_indices = sample_indices(validation_data, proportions_val)
    test_indices = sample_indices(test_data, proportions_test)

    # Apply the selected indices to create the new datasets
    training_data = (training_data[0][train_indices], (training_data[1][train_indices])-2)
    validation_data = (validation_data[0][val_indices], (validation_data[1][val_indices])-2)
    test_data = (test_data[0][test_indices], (test_data[1][test_indices])-2)
    """
    def func(pct, allvals):
      absolute = int(np.round(pct/100.*np.sum(allvals)))
      return f"{pct:.1f}%\n({absolute:d})"

    # DESCRIPTIVE ANALYSIS: oligomer states
    count = pd.Series(data[:,1]).value_counts()
    print(count)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Whole dataset")
    plt.show()
    count = pd.Series(training_data[1]).value_counts()
    print(count)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Training dataset")
    plt.show()
    count = pd.Series(test_data[1]).value_counts()
    print(count)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Testing dataset")
    plt.show()
    count = pd.Series(validation_data[1]).value_counts()
    print(count)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Validation dataset")
    plt.show()
    """
    training_data = (training_data[:,0], training_data[:,1])
    validation_data = (validation_data[:,0], validation_data[:,1])
    test_data = (test_data[:,0], test_data[:,1])
    """
    print(training_data[0].shape)
    print(training_data[1].shape)
    print(validation_data[0].shape)
    print(validation_data[1].shape)
    print(test_data[0].shape)
    print(test_data[1].shape)
    return (training_data, validation_data, test_data)

def vectorized_result(j):
    e = np.zeros((2, 1))
    e[j] = 1.0
    return e
"""
    e = np.zeros((2, 1))
    if j==3:
        e[1] = 1.0  # Label `2` as class `1`
    elif j==2:
        e[1] = 1.0
    elif j==4:
        e[1] = 1.0
    else:
        e[0] = 1.0  # All other digits (not 2) as class `0` (None)
    return e
"""
def load_data_wrapper():
    tr_d, va_d, te_d = load_data()
    # Convert string representations of lists into actual lists
    training_inputs = [eval(row) for row in tr_d[0]]
    # Convert to a NumPy array for easier manipulation
    training_inputs = np.array(training_inputs)
    training_inputs = [np.reshape(x, (28, 1)) for x in training_inputs]
    training_results = [vectorized_result(y) for y in tr_d[1]]
    training_data = list(zip(training_inputs, training_results))
    validation_inputs = [eval(row) for row in va_d[0]]
    validation_inputs = np.array(validation_inputs)
    validation_inputs = [np.reshape(x, (28, 1)) for x in validation_inputs]
    '''
    THEORY: ERROR IN ZIPPING PROCESS --> in python 3 list needs to be added to the zip command
    '''
    validation_data = list(zip(validation_inputs, va_d[1]))
    test_inputs = [eval(row) for row in te_d[0]]
    test_inputs = np.array(test_inputs)
    test_inputs = [np.reshape(x, (28, 1)) for x in test_inputs]
    test_data = list(zip(test_inputs, te_d[1]))
    return (training_data, validation_data, test_data)

training_data, validation_data, test_data = load_data_wrapper()

"""
modified_test_data = []
for x,y in test_data:
  if y == 3:
    modified_test_data.append((x,1))
  elif y == 2:
    modified_test_data.append((x,1))
  elif y == 4:
    modified_test_data.append((x,1))
  else:
    modified_test_data.append((x,0))
test_data = modified_test_data

modified_validation_data = []
for x,y in validation_data:
  if y == 3:
    modified_validation_data.append((x,1))
  elif y == 2:
    modified_validation_data.append((x,1))
  elif y == 4:
    modified_validation_data.append((x,2))
  else:
    modified_validation_data.append((x,0))
validation_data = modified_validation_data
"""
# Output of data generation and transformation:
print(f"Training data:")
print(training_data[0][0].shape)
print(training_data[0][1].shape)
print(f"\nValidation data:")
print(validation_data[0][0].shape)
print(f"\nTest data:\n")
print(test_data[0][0].shape)

# Neural Network applied to the coiled-coils data:
neurons = range(7,100+1)
neurons = [10]
epochs_total = []
accuracies_total = []
losses_total = []
max_accuracies = []
for neuron in neurons:
  NN = NeuralNetwork([28, neuron, 2])
  # Output for class definition:
  print(NN.layers)
  print(NN.sizes)
  print("Biases:")
  for bias in NN.biases:
      print(bias.shape)
  print("weights")
  for weigth in NN.weights:
      print(weigth.shape)
  # Output for feedforward architecture:
  epochs, accuracies, losses = NN.SGD(training_data, 30, 10, 3.0, test_data=test_data)
  epochs_total.append(epochs)
  accuracies_total.append(accuracies)
  losses_total.append(losses)
  max_accuracies.append(max(accuracies))

fig = plt.figure()
plt.plot(neurons, max_accuracies)
"""
neurons_summary = [i for i in range(7, 100+1, 5)]
max_accuracies_summary = []
print(neurons_summary)
for i in neurons_summary:
  print(f"Neuron: {i}")
  max_accuracies_summary.append(max_accuracies[neurons.index(i)])
print(max_accuracies_summary)
#plt.plot(neurons_summary, max_accuracies_summary, "red")
"""
plt.xlabel('Neurons in hidden layer')
plt.ylabel('Maximum accuracy')
plt.title('Accuracy vs Neurons in hidden layer')
plt.tight_layout()
plt.show()

# Plot the accuracy
fig1 = plt.figure()
for i in range(len(neurons)):
  plt.plot(range(epochs_total[i]), accuracies_total[i])
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.title('Accuracy vs Epochs')
#plt.legend(neurons, loc="center right",title="Neurons in hidden layer")
plt.tight_layout()
plt.show()

fig2 = plt.figure()
for i in range(len(neurons)):
  plt.plot(range(epochs_total[i]), losses_total[i])
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.title('Loss vs Epochs')
#plt.legend(neurons, loc="best",title="Neurons in hidden layer")
plt.tight_layout()
plt.show()

for i, ele in enumerate(training_data):
    x,y = ele
    y = np.argmax(y)
    training_data[i] = (x, y + 2)
