"""
Stacked RNN Architecture for Sequence Data Processing

Authors:
- Marieke van Vreeswijk
- Daniel Rademaker

Description:
This script contains the implementation of a stacked RNN architecture.
aimed at generating peptide sequence data.
"""

from rnn_module import RNNModule

class StackedRNN(nn.Module):
    """
    Implementation of a stacked RNN architecture.
    """

    def __init__(self, input_size, hidden_size, output_size, num_layers):
        """
        Initialize the stacked RNN architecture.
        
        :param input_size: Number of input features.
        :param hidden_size: Number of hidden layer features.
        :param output_size: Number of output features.
        :param num_layers: Number of RNN layers to stack.
        """
        super(StackedRNN, self).__init__()

        # If only one layer is specified, create a single RNN module with softmax activation
        if num_layers == 1:
            self.rnn_layers = [RNNModule(input_size, hidden_size, output_size, F.softmax)]
        else:
            self.rnn_layers = [RNNModule(input_size, hidden_size, hidden_size)]
            self.rnn_layers += [RNNModule(input_size, hidden_size, hidden_size) for _ in range(1, num_layers - 1)]
            self.rnn_layers += [RNNModule(hidden_size, hidden_size, output_size, F.softmax)]

    def forward(self, input_data):
        """
        Forward pass of the stacked RNN.
        
        :param input_data: Input tensor data.
        :return: Tensor representing the output of the stacked RNN.
        """
        for rnn_layer in self.rnn_layers:
            input_data = rnn_layer(input_data)
        return input_data

    def reset_memory(self):
        """
        Reset the internal memory of all RNN layers.
        """
        _ = [layer.reset_memory() for layer in self.rnn_layers]

    def eval(self):
        """
        Set all RNN layers to evaluation mode.
        """
        _ = [layer.eval() for layer in self.rnn_layers]

    def train(self):
        """
        Set all RNN layers to training mode.
        """
        _ = [layer.train() for layer in self.rnn_layers]

    def parameters(self):
        """
        Return parameters of all RNN layers.

        :return: List of all parameters from all RNN layers.
        """
        return [param for param in itertools.chain(*[layer.parameters() for layer in self.rnn_layers])]
