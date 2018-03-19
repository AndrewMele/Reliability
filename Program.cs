using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SoftwareReliability
{
    class Program
    {
        //Global Variables
        //static double[] component_reliabilities = { 0.99, 0.995, 0.98, 0.995};        //RAMS2014
        static double[] component_reliabilities = { 0.8, 0.75, 0.7, 0.72, 0.73 };   //Hoda's
        //static double[] component_reliabilities = { 0.99, 0.995, 0.98, 0.995, 0.8, 0.75, 0.7, 0.72, 0.73, 0.84 }; //10 State
        static int number = 5;
        static Matrix<double> correlations_matrix;
        static List<double> sigma = new List<double>();
        static List<int[]> CS = new List<int[]>();
            
        static void Main(string[] args)
        {
            //Macro to build Matrix
            var M = Matrix<double>.Build;
            //Double arrays based off of column_major notation which is used to create Matrices
            /*double[] correlations = {1.0000, 0.0094, -0.0023, -0.0020,              //RAMS2014
                                     0.0094, 1.0000,  0.0011,  0.0010,
                                    -0.0023, 0.0011,  1.0000,  0.0010,
                                    -0.0020, 0.0010,  0.0010,  1.0000}; */
            double[] correlations = {1.0000, -0.2800, 0.1000, -0.1500, 0.1500,    //Hoda's
                                      -0.2800,  1.0000, 0.2000,  0.3000,-0.1500,
                                       0.1000,  0.2000, 1.0000,  0.3200, 0.5000,
                                      -0.1500,  0.3000, 0.3200,  1.0000,-0.2200,
                                       0.1500, -0.1500, 0.5000, -0.2200, 1.0000};
            /*double[] correlations = {  1.0000,  0.2800, 0.1000, -0.0015, 0.1500, 0.0345, 0.0762, -0.0257, 0.0440, -0.0057,    //10 State
                                       0.2800,  1.0000, 0.2000,  0.3000,-0.0150, 0.0067,-0.0229,  0.0033, 0.1123,  0.0478,  
                                       0.1000,  0.2000, 1.0000,  0.3200, 0.0500,-0.0291, 0.0622,  0.0783, 0.0492, -0.0111,
                                      -0.0015,  0.3000, 0.3200,  1.0000,-0.0220, 0.0118, 0.0292,  0.0333, 0.0228,  0.0459,
                                       0.1500, -0.0150, 0.0500, -0.0220, 1.0000, 0.0234, 0.0456, -0.0654, 0.0122, -0.0555,
                                       0.0345,  0.0067,-0.0291,  0.0118, 0.0234, 1.0000, 0.0137, -0.0372, 0.0177,  0.0433,   
                                       0.0762, -0.0229, 0.0622,  0.0292, 0.0456, 0.0137, 1.0000,  0.0144, 0.0633, -0.0221,
                                      -0.0257,  0.0033, 0.0783,  0.0333,-0.0654,-0.0372, 0.0144,  1.0000, 0.0286,  0.0166,
                                       0.0440,  0.1123, 0.0492,  0.0228, 0.0122, 0.0177, 0.0633,  0.0286, 1.0000, -0.0223,
                                      -0.0057,  0.0478,-0.0111,  0.0459,-0.0555, 0.0433,-0.0221,  0.0166,-0.0223,  1.0000};*/

            //Creating Corresponding Matrices
            correlations_matrix = M.Dense(number, number, correlations);
            

            //10 State Cutset
            List<int[]> CS1 = new List<int[]>();
            CS1.Add(new int[] { 1, 6 });
            CS1.Add(new int[] { 4, 8 });
            CS1.Add(new int[] { 1, 5 });
            CS1.Add(new int[] { 4, 7 });
            CS1.Add(new int[] { 2, 5, 9 });
            CS1.Add(new int[] { 3, 5, 9 });
            CS1.Add(new int[] { 2, 8, 10 });
            CS1.Add(new int[] { 3, 6, 9 });
            CS1.Add(new int[] { 2, 7, 10 });
            CS1.Add(new int[] { 3, 8, 10 });
            CS1.Add(new int[] { 3, 7, 10 });
            CS1.Add(new int[] { 2, 6, 9 });
            CS1.Add(new int[] { 1, 7, 9, 10 });
            CS1.Add(new int[] { 4, 5, 9, 10 });
            CS1.Add(new int[] { 4, 6, 9, 10 });
            CS1.Add(new int[] { 1, 8, 9, 10 });
            //Hoda's Cutset
            List<int[]> CS2 = new List<int[]>();
            CS2.Add(new int[] {1,2});
            CS2.Add(new int[] {4,5});
            CS2.Add(new int[] {1,3,5});
            CS2.Add(new int[] {2,3,4});
            


            //Enumeration testing 
            number = 5;
            CS.Clear();
            CS.Add(new int[] {1,2});
            CS.Add(new int[] {4,5});
            CS.Add(new int[] {1,3,5});
            CS.Add(new int[] {2,3,4});
            List<int[]> Test = new List<int[]>();
            Test = RecursiveEnumerate(number);

            //Testing
            number = 5;
            //GenerateComponents(0.7, 0.99);
            //GenerateCorrelations(0.1);
            GenerateSigma();
            Console.Write("Component Reliabilities: < ");
            foreach (double d in component_reliabilities)
                Console.Write(d + " ");
            Console.WriteLine(">\n");
            Console.WriteLine("Correlations Matrix:");
            double[] correlations_array = correlations_matrix.AsColumnMajorArray();
            for (int i = 0; i < (number*number); i++)
            {
                if (i % number == 0)
                    Console.WriteLine("");
                Console.Write(correlations_array[i].ToString("#0.00000") + "  ");
            }
            Console.WriteLine("");

            FindReliability(number);
            double S5 = Simulation(1000000);
            double S6 = Simulation(10000000);
            S6 = 0;
            //Hold Results on screen
            Console.ReadLine();          
            
        }

        static double Simulation(int runs)
        {
            List<int> Path = new List<int>();
            Random Rand = new Random();
            int choice;
            int success = 0;

            for (int i = 0; i < runs; i++)
            {
                Path.Clear();
                for (int j = 0; j < number; j++)
                {
                    choice = Rand.Next(0,2);
                    Path.Add(choice);
                    if (choice == 0 && ContainsCutset(Path.ToArray()))
                        break;
                }
                if (Path.Count == number)
                    success += 1;

            }
            return (double) success / (double) runs; 
        }

        static void PlotConvergence()
        {
            //success_percentages[0] correlates to the first set of data from m = 0 -> 1.0
            List<double[]> success_percentages = new List<double[]>();
            for (number = 10; number <= 10; number++)
                success_percentages.Add(ConvergenceTest(0.1,1000));
            Console.WriteLine("All Percentages Added");


        }
        static double[] ConvergenceTest(double step_size, int runs)
        {
            List<double> percentages = new List<double>(); //List to hold percentages of all steps from 0 -> 1.0
            double result;      //resulting percentage after x amount of runs
            double sum;           //Total Reliability of the system, if this value exceeds 1 or is lower than 0, we know there was an error
            bool neg_flag;      //Flag that determines whether or not a set of probabilities had a negative value contained within it
            int success_counter;//Counts the number of successful operations
            
            //Enumerate Final States
            List<int[]> Final_States = RecursiveEnumerate(number);

            for (double m = 0; m < 1.0; m+=step_size)   //Loop to increment multiplier after x amount of runs
            {
                //reset counters after step size increment 
                success_counter = 0;

                for (int x = 0; x < runs; x++)          //Loop to run through operation x amount of times
                {
                    GenerateComponents(0.9,1.0);
                    GenerateCorrelations(m);
                    GenerateSigma();
                    sum = 0;
                    neg_flag = false;
                    foreach (int[] i in Final_States)
                    {
                        sum += TotalProbability(i);
                        if (TotalProbability(i) < 0)
                            neg_flag = true;
                    }
                    if (!neg_flag && Math.Abs(1.0-sum) <= 0.00001)
                        success_counter++;

                }
                result = (double)success_counter/(double)runs;
                percentages.Add(result);
                Console.WriteLine("number = {0}\t m = {1}\t{2}\t complete...", number.ToString(), m.ToString("#0.00"), (result).ToString("###.00%"));  
            }
            return percentages.ToArray();
        }
        struct CorrelationPairs
        {
            public int x;
            public int y;
            public double value;

            public CorrelationPairs(int xx, int yy, double v)
            {
                x = xx;
                y = yy;
                value = v;
            }
        };

        //Functions

        //Function which randomly creates a number of components based on the min and max reliability given
        static void GenerateComponents(double min, double max)
        {
            List<double> temp = new List<double>();
            Random R = new Random();
            for (int i = 0; i < number; i++)
                temp.Add(R.NextDouble() * (max-min) + min);
            component_reliabilities = temp.ToArray();
        }

        static void GenerateCorrelations(double multi = -1)
        {
            List<CorrelationPairs> bounds = CorrelationBounds(multi);
            List<double> temp = new List<double>();
            for (int i = 0; i < number; i++)
            {
                for (int j = 0; j < number; j++)
                {
                    if (i == j)
                        temp.Add(1);
                    else if (j < i)
                        temp.Add(bounds.Find(z => z.x == j && z.y == i).value);
                    else
                        temp.Add(bounds.Find(z => z.x == i && z.y == j).value);
                }
            }
            
           correlations_matrix = Matrix<double>.Build.Dense(number,number,temp.ToArray());

        }

        static void GenerateSigma()
        {
            sigma.Clear();
            foreach (double d in component_reliabilities)
            {
                sigma.Add(Math.Sqrt(d * (1 - d)));
            };
        }

        //Function which finds the maximum and minimum bounds for correlation given each combination
        static List<CorrelationPairs> CorrelationBounds(double multi = -1)
        {
            List<CorrelationPairs> temp = new List<CorrelationPairs>();
            Random R = new Random();
            double min = 0;
            double max;
            double bound;
            double sqrt1;
            double sqrt2;
            //Find Min Bound Correlations
            for (int i = 0; i < component_reliabilities.Count() - 1; i++)
            {
                for (int j = i + 1; j < component_reliabilities.Count(); j++)
                {
                    //I think these equations are wrong
                    sqrt1 = -Math.Sqrt(component_reliabilities[i] * component_reliabilities[j] / ((1 - component_reliabilities[i]) * (1 - component_reliabilities[j])));
                    sqrt2 = -Math.Sqrt(((1 - component_reliabilities[i]) * (1 - component_reliabilities[j])) / (component_reliabilities[i] * component_reliabilities[j]));
                    min = Math.Max(sqrt1, sqrt2);
                    if (multi != -1)
                        min *= multi;
                    sqrt1 = Math.Sqrt((component_reliabilities[i] * (1 - component_reliabilities[j])) / (component_reliabilities[j] * (1 - component_reliabilities[i])));
                    sqrt2 = Math.Sqrt((component_reliabilities[j] * (1 - component_reliabilities[i])) / (component_reliabilities[i] * (1 - component_reliabilities[j])));
                    max = Math.Min(sqrt1, sqrt2);
                    if (multi != -1)
                        max *= multi;
                    //Overwrite equations with hard coded values
                    //min = -0.01;
                    //max = 0.1;
                    bound = R.NextDouble() * (max - min) + min;
                    /* if (multi != -1)
                        Console.WriteLine("Multiplier = {0}", multi.ToString("#0.00"));
                    Console.WriteLine("u[{0}]u[{1}] =\t {2} < x < {3}\tx= {4}", i, j, min.ToString("#0.0000000000"), max.ToString("#0.0000000000"), bound.ToString("#0.0000000000"));
                    */temp.Add(new CorrelationPairs(i, j, bound));
                }
            }
            return temp;
        }

        //Function that should reproduce the results of the RAMS2014 Paper (Make sure correlations matrix and component reliabilities are changed)
        static void RAMS2014()
        {
            List<int[]> All_States = RecursiveEnumerate(4);
            foreach (int[] i in All_States)
            {
                Console.Write("Pr(");
                foreach (int ii in i)
                    Console.Write(ii);
                Console.WriteLine(") = " + TotalProbability(i).ToString("#0.############"));
            }
            Console.WriteLine("\nE{R}\t = " + (TotalProbability(new int[] { 1, 1, 1, 1 }) + TotalProbability(new int[] { 1, 1, 1, 0 }) + TotalProbability(new int[] { 1, 1, 0, 1 }) + TotalProbability(new int[] { 1, 1, 0, 0 }) + TotalProbability(new int[] { 1, 0, 1, 1 }) + TotalProbability(new int[] { 0, 1, 1, 1 })).ToString("#0.############"));
            Console.WriteLine("E{FS}\t = " + (TotalProbability(new int[] { 1, 0, 1, 0 }) + TotalProbability(new int[] { 1, 0, 0, 1 }) + TotalProbability(new int[] { 0, 1, 1, 0 }) + TotalProbability(new int[] { 0, 1, 0, 1 })).ToString("#0.############"));
            Console.WriteLine("E{UF}\t = " + (TotalProbability(new int[] { 1, 0, 0, 0 }) + TotalProbability(new int[] { 0, 1, 0, 0 }) + TotalProbability(new int[] { 0, 0, 1, 1 }) + TotalProbability(new int[] { 0, 0, 1, 0 }) + TotalProbability(new int[] { 0, 0, 0, 1 }) + TotalProbability(new int[] { 0, 0, 0, 0 })).ToString("#0.############"));
        }

        //Function which calculates bridge reliability based on n number of components and cutsets CS
        static double FindReliability(int n)
        {
            //Enumerate All Successful States
            List<int[]> Paths = RecursiveEnumerate(n);
            double R = 0;
            //Loop through and see which Leaf Nodes are success
            foreach (int[] i in Paths)
            {
                    Console.Write("Pr(");
                    foreach (int ii in i)
                        Console.Write(ii);
                    Console.WriteLine(") \t= " + TotalProbability(i).ToString("#0.0000000000"));
                    R += TotalProbability(i);
                    Console.WriteLine("Reliability \t= " + R.ToString("#0.0000000000") + "\n");
            }
            return R;
        }

        //Recursive Dynamic Approach to Enumerating all possible states, defaults to only adding Full Path and not partial
        static List<int[]> RecursiveEnumerate(int n)
        {
            int[] counters = new int[n];
            for (int i = 0; i < n; i++)
                counters[i] = 0;
            List<int[]> Paths = new List<int[]>();

            Paths = RecursiveEnumerate(0, counters, Paths);
            return Paths;          
        }

        //Second Form of Recursive Function
        static List<int[]> RecursiveEnumerate(int level, int[] counters, List<int[]> Paths)
        {
            List<int> temp_state = new List<int>();
            
            for (int i = 0; i < 2; i++)
            {
                temp_state.Clear();
                counters[level] = i;
                for (int a = 0; a <= level; a++)
                    temp_state.Add(counters[a]);
                if (level != counters.Length - 1 && !ContainsCutset(temp_state.ToArray()))
                    Paths = RecursiveEnumerate(level + 1, counters, Paths);
                if (level == counters.Length - 1 && !ContainsCutset(temp_state.ToArray()))
                    Paths.Add(temp_state.ToArray());
            }
            return Paths;
        }

        //Check to see if the integer array is a cutset, returns true if the query is a cutset, returns false otherwise
        static bool ContainsCutset(int[] query)
        {
            bool safe;
            foreach (int[] cut in CS)            //Loop through each cutset and compare to query
            {
                safe = false;
                if (query.Count() >= cut.Max())     //See if check is applicable
                {
                    foreach (int i in cut)
                    {
                        if (query[i - 1] != 0)
                            safe = true;
                    }
                    if (!safe)
                        return true;
                }
                    
            }
            return false;
        }

        //Static Enumeration of all States with n=5 Final option only adds arrays of count n (Looking to make this dynamic in the future... Possibly Cartesian Product?)
        static List<int[]> EnumerateAllStates(bool Final = true)
        {
            List<int> temp_state = new List<int>();
            List<int[]> All_States = new List<int[]>();
            List<int[]> Final_States = new List<int[]>();

            //Enumerate all states
            for (int a = 0; a < 2; a++)
            {
                temp_state.Add(a);
                All_States.Add(temp_state.ToArray());
                temp_state.Clear();
                for (int b = 0; b < 2; b++)
                {
                    temp_state.Add(a);
                    temp_state.Add(b);
                    All_States.Add(temp_state.ToArray());
                    temp_state.Clear();
                    for (int c = 0; c < 2; c++)
                    {
                        temp_state.Add(a);
                        temp_state.Add(b);
                        temp_state.Add(c);
                        All_States.Add(temp_state.ToArray());
                        temp_state.Clear();
                        for (int d = 0; d < 2; d++)
                        {
                            temp_state.Add(a);
                            temp_state.Add(b);
                            temp_state.Add(c);
                            temp_state.Add(d);
                            All_States.Add(temp_state.ToArray());
                            temp_state.Clear();
                            for (int e = 0; e < 2; e++)
                            {
                                temp_state.Add(a);
                                temp_state.Add(b);
                                temp_state.Add(c);
                                temp_state.Add(d);
                                temp_state.Add(e);
                                All_States.Add(temp_state.ToArray());
                                Final_States.Add(temp_state.ToArray());
                                temp_state.Clear();
                            }
                        }
                    }
                }
            }
            if (Final)
                return Final_States;
            else
                return All_States;
        }

        //Random Tree Traversal with Cutsets (n = max level of tree)
        static void RTwithCS(int n)
        {
            List<int> path = new List<int>();
            Random r = new Random();
            int choice;
            Console.WriteLine("Starting Random Traversal with Cutsets");
            while (true)
            {
                choice = r.Next(0, 2);
                Console.Write("Choose " + choice + "\tPath: ");
                path.Add(choice);
                for (int i=0; i<path.Count; i++)
                {
                    if (i == path.Count - 1)
                        Console.WriteLine(path[i]);
                    else
                        Console.Write(path[i] + ",");
                }

                if (ContainsCutset(path.ToArray()))
                {
                    Console.WriteLine("Cutset Determined!\nAborting Traversal...\n");
                    break;
                }
                if (path.Count == n && !ContainsCutset(path.ToArray()))
                {
                    Console.WriteLine("Success!\n");
                    break;
                }
            }
            Console.ReadLine();
        }

        //Simple Random Traversal checking k_out_of_n standards, cons value determines if the failures must be consecutive
        static void RandomTraversal(int k, int n, bool cons = true)
        {
            List<int> temp = new List<int>();
            Random r = new Random();
            int choice;
            int count;
            int conscount;
            Console.WriteLine("Starting Random Traversal");
            while (true)
            {
                count = 0;
                conscount = 1;
                choice = r.Next(0, 2);
                Console.WriteLine("Choose " + choice);
                temp.Add(choice);
                Console.Write("Path: {");
                foreach (int i in temp)
                {
                    Console.Write(i+",");
                    if (i == 0)
                        count++;
                }
                if (cons)
                {
                    for (int i = 1; i < temp.Count; i++)
                    {
                        if (temp[i] == 0 && temp[i - 1] == 0)
                            conscount++;
                        else
                            conscount = 1;
                    }
                }
                Console.WriteLine("}");
                if (cons)
                {
                    if (conscount >= k || temp.Count == n)
                        break;
                }

                else
                {
                    if (count >= k || temp.Count == n)
                        break;
                }
                
            }
            if (cons && conscount >= k)
                Console.WriteLine("Consecutive Failure Threshold Reached!\nEnding random traversal...\n");
            else if (!cons && count >= k)
                Console.WriteLine("Failure Threshold Reached!\nEnding random traversal...\n");
            else
                Console.WriteLine("Success!");
            Console.ReadLine();

        }


        //Prints out all acceptable or not acceptable states using k_out_of_n function, OK value determines which set is returned
        static List<int[]> k_out_of_n(int k, int n, bool OK = true)
        {
            List<int[]> All_States = new List<int[]>();
            List<int[]> OK_States = new List<int[]>();
            List<int[]> BAD_States = new List<int[]>();
            List<int> temp_state = new List<int>();
            int count;

            All_States = RecursiveEnumerate(n);
            
            foreach (int[] i in All_States)
            {
                count = 0;
                foreach (int j in i)
                    if (j == 0)
                        count++;
                if (count < k)
                    OK_States.Add(i);
                else
                    BAD_States.Add(i);
            }

            foreach (int[] i in OK_States)
            { 
                foreach (int j in i)
                    Console.Write(j);
                Console.WriteLine();
            }
            Console.ReadLine();

            if (OK)
                return OK_States;
            else
                return BAD_States;
            
        }
       
        //Prints out Total Probability for all 32 Combinations
        static void PrintTotal32()
        {
            //Print out Component Reliabilites and Correlations Matrix
            Console.Write("Component Reliabilities: {");
            component_reliabilities.ToList().ForEach(i => Console.Write(i.ToString() + ", "));
            Console.WriteLine("}\n");
            Console.WriteLine("Correlations Matrix:");
            Console.WriteLine(correlations_matrix.ToString());
            List<int[]> FinalStates = RecursiveEnumerate(5);
            foreach (int[] i in FinalStates)
            {
                Console.Write("Pr(");
                foreach (int j in i)
                    Console.Write(j);
                Console.WriteLine(")\t=\t" + TotalProbability(i).ToString("#0.############"));
            }
         }

        //TotalProbability takes in an integer array which corresponds to the binary string of successes and failures and returns the appropriate probability with deconditoning
        //Example: 1010 = Fourth State being failure given first state was success and second state was failure and third state was success
        static double TotalProbability(int[] T1, bool print = false)
        {
            //Create List to hold the intermediate Terms and temp variable to create terms as well as final term variable and Result variable
            List<int[]> T = new List<int[]>();
            int[] temp;
            int FinalTerm = T1.First();
            double R = 1;

            //Add the rest of the terms
            for (int i = 0; i < T1.Length - 1; i++)
            {
                temp = T1;
                Array.Resize(ref temp, temp.Length - i);
                T.Add(temp);
            }

            //Multiply each ConditionalProbability in T with R
            foreach (int[] i in T)
                R *= ConditionalProbability(i,print);

            //Multiply FinalTerm
            if (FinalTerm == 1)
                R *= component_reliabilities[0];
            else
                R *= (1 - component_reliabilities[0]);

            return R;
        }

        //ConditionalProbability takes in an integer array which corresponds to the binary string of successes and failures and returns the appropriate probability
        //Example: 1010 = Fourth State being failure given first state was success and second state was failure and third state was success
        static double ConditionalProbability(int[] B, bool print = false)
        {
            double R = 0;

            if (B.Last() == 0 && print)
                Console.Write("1 - (");

            R = component_reliabilities[B.Length - 1];

            if (print)
                Console.Write("u[" + B.Length + "] " + component_reliabilities[B.Length - 1] + " + ");

            for (int i = 0; i < B.Length - 1; i++)
            {
                if (B[i] == 1)
                {
                    R += (1 - component_reliabilities[i]) * b(B.Length, i + 1);
                    if (print)
                        Console.Write("(1 - u[" + (i+1) + "])b(" + (B.Length) + ", " + (i + 1) + ")");
                }
                else
                {
                    R += (0 - component_reliabilities[i]) * b(B.Length, i + 1);
                    if (print)
                        Console.Write("(0 - u[" + (i + 1) + "])b(" + (B.Length) + ", " + (i + 1) + ")");
                }
                if (print)
                {
                    if (i != B.Length - 2)
                      Console.Write(" + ");
                    else
                      Console.Write(") = ");
                }
            }
            if (print)
            {
                if (B.Last() == 1)
                    Console.WriteLine(R);
                else
                    Console.WriteLine(1 - R);
            }
            
            if (B.Last() == 1)
                return R;
            else
                return 1 - R;
        }

        //The b function used in the ConditionalProbability function
        static double b(int i, int j, bool print = false)
        {
            //Create temporary list to make Matrices out of
            List<double> temp = new List<double>();
            for (int k = 0; k < i-1; k++)
            {
                temp.AddRange(correlations_matrix.Column(k).Take(i-1));
            }
            //Create First Matrix out of temp array
            Matrix<double> M1 = Matrix<double>.Build.Dense(i - 1, i - 1, temp.ToArray());
            if (print)
                Console.WriteLine("M1 Matrix Before Sigma:\n" + M1.ToString());
            //Multiply First Matrix values by the appropriate sigma values
            for (int k = 0; k < i - 1; k++)
                for (int h = 0; h < i - 1; h++)
                    M1[k, h] = M1[k, h] * sigma[k] * sigma[h];
            //Create Second Matrix (must clear temp list)
            if (print)
                Console.WriteLine("M1 Matrix After Sigma:\n" + M1.ToString());
            temp.Clear();
            temp.AddRange(correlations_matrix.Column(i - 1).Take(i - 1));
            Matrix<double> M2;
            M2 = Matrix<double>.Build.Dense(i - 1, 1, temp.ToArray());
            if (print)
                Console.WriteLine("M2 Matrix Before Sigma:\n" + M2.ToString());
            //Multiply Second Matrix by appropriate sigma values
            for (int k = 0; k < i - 1; k++)
                    M2[k, 0] = M2[k, 0] * sigma[k] * sigma[i - 1];
            //Put Results into new Matrix and return the value
            if (print)
            {
                Console.WriteLine("M2 Matrix After Sigma:\n" + M2.ToString());
                Console.WriteLine("M1 Inverse Matrix:\n" + M1.Inverse().ToString());
            }
            Matrix<double> R = M1.Inverse() * M2;
            if (print)
                Console.WriteLine("M1 Inverse Matrix * M2 Matrix:\n" + R.ToString());
            return R[j-1,0];
            
        }
    }
}
