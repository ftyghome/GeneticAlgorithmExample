## Sample of the Genetic Algorithm
The program solves a simple problem, that is minimizing the value of a simple function $f(x1,x2)=21.5+x_1*sin(4\pi x_1)+x_2*sin(20\pi x_2)$

![f(x1,x2)](README.assets/f(x1,x2)-1616050406355.png)

### Usage

Modify the main function allows you to change the parameters of the algorithm

The constructor of class GeneticAlgo accepts seven parameters:



| Name                | Explanation                                                  |
| ------------------- | ------------------------------------------------------------ |
| Func*               | The target function the program want to get its maximum      |
| BitConvertTools& _a | Related to the domain of definition and the precision of the $x_1$ |
| BitConvertTools& _b | Related to the domain of definition and the precision of the $x_2$ |
| int _N              | The population of ants                                       |
| double _pc          | The possibility of crossover                                 |
| double _pm          | The possibility of metamorphosis                             |
| int _Gmax           | Number of iterations                                         |

