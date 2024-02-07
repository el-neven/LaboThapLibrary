from Port.Mass_connector import Mass_connector

class Cycle:
    def __init__(self, number_of_points, list_of_components):
        """
        Initialize a Cycle object.

        Args:
            number_of_points (int): The number of states on the P&ID, starting from 0.
            list_of_components (list): A list of component objects.

        Attributes:
            defined (bool): Indicates whether the cycle is defined or not.
            n (int): The number of points.
            compo_list (list): The list of components.
            point (list): The list of mass connectors.

        """
        self.defined = False 
        self.n = number_of_points
        self.compo_list = list_of_components
        
        self.point = []
        for i in range(self.n):
            self.point.append(Mass_connector())

        self.check_parametrized() #Check if the parameters of the components are known
        

    def check_parametrized(self):
        self.defined = True # Assumption
        
        for i in range(len(self.compo_list)):
            if self.compo_list[i].parametrized:
                pass
            else:
                self.defined = False
                print(f"Error: Parameters of {self.compo_list[i]} not completely known")
    
    def solve(self):
        if self.defined:
            n_solved = 0
            n_it = 0
            while n_solved<len(self.compo_list) and n_it<len(self.compo_list)+1:
                for i in range (len(self.compo_list)):
                    # Search the next component calculable
                    if not self.compo_list[i].calculable:
                        self.compo_list[i].check_calculable()

                    if self.compo_list[i].calculable and not self.compo_list[i].defined: #If the component is calculable and not defined then we solve it
                        self.compo_list[i].solve()
                        n_solved +=1

                    else:
                        pass
                
                n_it +=1
            if n_solved == len(self.compo_list):
                print("Cycle solved")
                
            else:
                print("Error: Cycle not correctly defined")

        else:
            print("Error: Parameters of the component not completely known")
            
        
        
        
    
