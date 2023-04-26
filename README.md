# Shors-Algorithm
An implementation of Shor's quantum factoring algorithm on the number 15:

 - Uses IBM's qiskit python API for quantum circuit contruction and simulation.
qiskit can be found here: https://qiskit.org/

 - Successfully finds the factors of 15 (3 and 5) using the properties of quantum superposition and phase estimation.

Included are a .py file containing the raw code, a Jupyter Notebook file (.ipynb) which contains both the code and the output,
and a slide deck which contains an explaination of the algorithm, its implentation, and the underlying math.


Examples of different amod15 gate circuits:

![amod15 gate circuits.](https://github.com/aklein4/Shors-Algorithm/blob/main/Example%20Images/amod15_gates.jpg)

Example of a QFT inversion gate circuit with 8 qubits:

![QFT Inversion Circuit.](https://github.com/aklein4/Shors-Algorithm/blob/main/Example%20Images/n8_QFT_inversion.jpg)

The full circuit outline:

![Diagram of the entire cirtuit outline.](https://github.com/aklein4/Shors-Algorithm/blob/main/Example%20Images/full_circuit.jpg)

Example of qubit readings after QFT Inversion:

![Graph showing probabilities of qubit readings after inversion.](https://github.com/aklein4/Shors-Algorithm/blob/main/Example%20Images/QFT_inversion_output.jpg)
