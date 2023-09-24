import torch
import itertools
import torch.nn as nn
import torch.nn.functional as F

class RNNModule(nn.Module):
    """
    Implementation of a single RNN cell module.
    """

    def __init__(self, input_size, hidden_size, output_size, output_activation=F.relu):
        """
        Initialize the RNN module.
        
        :param input_size: Number of input features.
        :param hidden_size: Number of hidden layer features.
        :param output_size: Number of output features.
        :param output_activation: Activation function to apply on the output.
        """
        super(RNNModule, self).__init__()

        concatenated_size = input_size + hidden_size
        self.memory_size = hidden_size

        self.input_to_hidden = nn.Linear(concatenated_size, hidden_size)
        self.hidden_to_hidden = nn.Linear(hidden_size, hidden_size)
        self.hidden_to_output = nn.Linear(hidden_size, output_size)
        self.hidden_to_memory = nn.Linear(hidden_size, hidden_size)
        self.dropout = nn.Dropout(p=0.1)
        self.output_activation = output_activation

        self.reset_memory() 

    def reset_memory(self):
        """
        Reset the internal memory of the RNN module.
        """
        self.memory = torch.zeros(self.memory_size)

    def forward(self, input_data):
        """
        Forward pass of the RNN module.
        
        :param input_data: Input tensor data.
        :return: Tensor representing the output of the module.
        """
        concatenated_data = torch.cat([input_data, torch.tanh(self.memory)], 0)
        hidden_data = F.relu(self.input_to_hidden(concatenated_data))
        hidden_data = F.relu(self.hidden_to_hidden(hidden_data))
        hidden_data = self.dropout(hidden_data)
        output_data = self.hidden_to_output(hidden_data)
        self.memory += self.hidden_to_memory(hidden_data)

        return self.output_activation(output_data)
