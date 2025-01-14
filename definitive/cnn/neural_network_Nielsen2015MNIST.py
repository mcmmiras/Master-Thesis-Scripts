# -*- coding: utf-8 -*-
"""
This script has been adapted using the original NN script developed by 
Nielsen in 2015.

REFERENCES:
- http://neuralnetworksanddeeplearning.com/chap1.html
- https://github.com/nchen909/NN-by-hand
"""


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

# To calculate the relu function of the NN
def relu(x):
    return np.maximum(0, x)

# To calculate the derivative of the relu function, used for backpropagation purposes
def relu_prime(x):
    return np.where(x > 0, 1, 0)

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
        # Note that the variable l in the loop below is used a little
        # differently to the notation in Chapter 2 of the book.  Here,
        # l = 1 means the last layer of neurons, l = 2 is the
        # second-last layer, and so on.  It's a renumbering of the
        # scheme in the book, used here to take advantage of the fact
        # that Python can use negative indices in lists.
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

# LOADING OF THE MNIST DATASET

# Libraries
import pickle
import gzip
import numpy as np

def load_data():
    f = gzip.open('/content/drive/MyDrive/Semestre 3 (2024-2025)/lab_shared/mnist.pkl.gz', 'rb')
    training_data, validation_data, test_data = pickle.load(f,encoding="latin1")
    """
    # Randomly shuffle the indices for each dataset
    train_indices = np.random.choice(training_data[0].shape[0], size=50000, replace=False)
    val_indices = np.random.choice(validation_data[0].shape[0], size=10000, replace=False)
    test_indices = np.random.choice(test_data[0].shape[0], size=10000, replace=False)

    # Apply the indices to randomly select the data
    training_data = (training_data[0][train_indices], training_data[1][train_indices])
    validation_data = (validation_data[0][val_indices], validation_data[1][val_indices])
    test_data = (test_data[0][test_indices], test_data[1][test_indices])
    f.close()
    """
    proportions_train = {
        0: 0,
        1: 0,
        2: 239/600,
        3: 126/600,
        4: 184/600,
        5: 13/600,
        6: 32/600,
        7: 4/600,
        8: 2/600,
        9: 0
    }
    proportions_test = {
        0: 0,
        1: 0,
        2: 48/129,
        3: 30/129,
        4: 39/129,
        5: 4/129,
        6: 6/129,
        7: 2/129,
        8: 0,
        9: 0
    }
    proportions_val = {
        0: 0,
        1: 0,
        2: 43/129,
        3: 32/129,
        4: 42/129,
        5: 2/129,
        6: 5/129,
        7: 5/129,
        8: 0,
        9: 0
    }
    """
    # Randomly shuffle the indices for each dataset
    train_indices = np.random.choice(training_data[0].shape[0], size=600, replace=False)
    val_indices = np.random.choice(validation_data[0].shape[0], size=129, replace=False)
    test_indices = np.random.choice(test_data[0].shape[0], size=129, replace=False)

    # Apply the indices to randomly select the data
    training_data = (training_data[0][train_indices], training_data[1][train_indices])
    validation_data = (validation_data[0][val_indices], validation_data[1][val_indices])
    test_data = (test_data[0][test_indices], test_data[1][test_indices])
    f.close()
    """
    # Function to sample indices based on class proportions
    def sample_indices(data, proportions, num_samples):
      x, y = data  # x = images, y = labels
      selected_indices = []

      # Calculate the total number of samples we need to pick for each class
      total_samples = num_samples

      for class_label, proportion in proportions.items():
          # Get the indices for the current class
          class_indices = np.where(y == class_label)[0]

          # Calculate how many samples we need to pick from this class
          num_class_samples = int(np.round(proportion * total_samples))
          num_class_samples = int(proportion*total_samples)
          if len(class_indices) > num_class_samples:
              # Randomly select the required number of samples for this class
              selected_class_indices = np.random.choice(class_indices, num_class_samples, replace=False)
          else:
              # If the class doesn't have enough samples, take all of them
              selected_class_indices = class_indices

          # Add the selected class indices to the list
          selected_indices.extend(selected_class_indices)

      return np.array(selected_indices).astype(int)

    # Apply sampling based on proportions
    train_indices = sample_indices(training_data, proportions_train, 600)
    test_indices = sample_indices(test_data, proportions_test, 129)
    val_indices = sample_indices(validation_data, proportions_val, 129)

    # Apply the indices to create the new datasets
    training_data = (training_data[0][train_indices], (training_data[1][train_indices])-2)
    test_data = (test_data[0][test_indices], (test_data[1][test_indices])-2)
    validation_data = (validation_data[0][val_indices], (validation_data[1][val_indices])-2)
    f.close()

    # DESCRIPTIVE ANALYSIS: digits
    def func(pct, allvals):
      absolute = int(np.round(pct/100.*np.sum(allvals)))
      return f"{pct:.1f}%\n({absolute:d})"

    count = pd.Series(training_data[1]).value_counts() + pd.Series(validation_data[1]).value_counts() + pd.Series(test_data[1]).value_counts()
    count = count.fillna(0)
    print(count)
    plt.pie(count, autopct="%1.1f%%")
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Whole dataset")
    plt.show()

    count = pd.Series(training_data[1]).value_counts() + pd.Series(validation_data[1]).value_counts() + pd.Series(test_data[1]).value_counts()
    count = count.fillna(0)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Whole dataset")
    plt.show()

    count = pd.Series(training_data[1]).value_counts()
    count = count.fillna(0)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Training dataset")
    plt.show()

    count = pd.Series(test_data[1]).value_counts()
    count = count.fillna(0)
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Testing dataset")
    plt.show()
    count = pd.Series(validation_data[1]).value_counts()
    plt.pie(count, autopct=lambda pct: func(pct, count))
    plt.legend(labels=count.index, loc="upper right")
    plt.title("Validation dataset")
    plt.show()

    print(training_data[0].shape)
    print(training_data[1].shape)
    print(validation_data[0].shape)
    print(validation_data[1].shape)
    print(test_data[0].shape)
    print(test_data[1].shape)
    return (training_data, validation_data, test_data)

def load_data_wrapper():
    tr_d, va_d, te_d = load_data()
    training_inputs = [np.reshape(x, (784, 1)) for x in tr_d[0]]
    training_results = [vectorized_result(y) for y in tr_d[1]]
    training_data = list(zip(training_inputs, training_results))

    validation_inputs = [np.reshape(x, (784, 1)) for x in va_d[0]]
    validation_data = list(zip(validation_inputs, va_d[1]))
    test_inputs = [np.reshape(x, (784, 1)) for x in te_d[0]]
    test_data = list(zip(test_inputs, te_d[1]))
    return (training_data, validation_data, test_data)

def vectorized_result(j):
    e = np.zeros((7, 1))
    e[j] = 1.0
    return e

# EXECUTION
training_data, validation_data, test_data = load_data_wrapper()

# Output of data generation and transformation:
print(f"Training data:")
print(training_data[0][0].shape)
print(training_data[0][1].shape)
print(f"\nValidation data:")
print(validation_data[0][0].shape)
print(validation_data[0][1].shape)
print(f"\nTest data:\n")
print(test_data[0][0].shape)
print(test_data[0][1].shape)

# Neural network:
neurons= range(7,100+1)
neurons = [10]
epochs_total = []
accuracies_total = []
losses_total = []
max_accuracies = []
for neuron in neurons:
  NN = NeuralNetwork([784, neuron, 7])
  # Output for class definition:
  print(NN.layers)
  print(NN.sizes)
  print("Biases:")
  for bias in NN.biases:
      print(bias.shape)
  print("weights")
  for weigth in NN.weights:
      print(weigth.shape)
  # Output for the feedforward architecture of the NN:
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
plt.plot(neurons_summary, max_accuracies_summary, "red")
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
