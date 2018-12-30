using System;
using System.Collections;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SoftwareReliability
{
    class Program
    {
        //Global Variables
        static int number;
        static double[] component_reliabilities;
        static Matrix<double> correlations_matrix;
        static List<double> sigma = new List<double>();
        static List<int[]> CS = new List<int[]>();
        static List<int[]> PS = new List<int[]>();
        static List<double[]> Bm = new List<double[]>();
        //static double[] component_reliabilities = {0.95, 0.8, 0.93};                    //2 out of 3
        //static double[] component_reliabilities = { 0.99, 0.995, 0.98, 0.995};        //RAMS2014
        //component_reliabilities = { 0.8, 0.75, 0.7, 0.72, 0.73 };   //Hoda's
        //static double[] component_reliabilities = { 0.99, 0.995, 0.98, 0.995, 0.8, 0.75, 0.7, 0.72, 0.73, 0.84 }; //10 State

        static List<PTrace> ConditionalPErrors = new List<PTrace>();
        static List<PTrace> TotalPErrors = new List<PTrace>();

        static void Main(string[] args)
        {
            //Setup stopwatch
            DateTime startTime, endTime;
            double elpasedMilliseconds;
            startTime = DateTime.Now;

            //Build environment
            Build("8p");
            

            int i = 1;
            double S;
            if (i == 0)
            {    
                FindReliability();
            }
            else
            {
                S = Simulation(1000000);
                Console.WriteLine("Simulation Estimate with 1,000,000 runs: " + S.ToString("#0.0000000000"));
            }

            //Hold Results on screen
            endTime = DateTime.Now;
            elpasedMilliseconds = ((TimeSpan)(endTime - startTime)).TotalMilliseconds;
            Console.WriteLine("Task took " + elpasedMilliseconds.ToString("#0.00") + " ms.");
            //ExportBuild("22bn");         
        }
        
        //Exports current enviornment to a text file which includes: component reliabilities and correlations
        static void ExportBuild(string file_name)
        {
            string path = @"Z:\" + file_name + ".txt";
            string components = "{";
            string correlations_string = "{";
            double [] correlations = correlations_matrix.ToColumnMajorArray();
            //Start File
            File.AppendAllLines(path, new [] {file_name + "Build", "Component Reliabilities:"});
            //Loop through component reliabilities
            for (int i = 0; i < component_reliabilities.Count(); i++)
            {
                if (i == component_reliabilities.Count() - 1)
                    components += component_reliabilities[i].ToString("#0.00000") + " }";
                else
                    components += component_reliabilities[i].ToString("#0.00000") + ", ";
            }
            //Append components
            File.AppendAllLines(path, new [] {components, "Correlations:"});
            //Loop through correlations
            for (int i = 0; i < correlations.Count(); i++)
            {
                if (i == correlations.Count() - 1)
                    correlations_string += correlations[i].ToString("#0.00000") + " }";
                else if (i % number == 0)
                    correlations_string += correlations[i].ToString("#0.00000") + ",\n";
                else
                    correlations_string += correlations[i].ToString("#0.00000") + ", ";
            }
            //Append correlations
            File.AppendAllLines(path, new [] {correlations_string});
        }

        //Function which uses a Random Number Generator and a data set to sample from the data set and generate a double within the bounds of the sample
        static double PseudoSample(double[] data)
        {
            Random r = new Random();        //Random Number Generator to choose random correlations from the data set
            int a;                          //Temp variable to hold random number between 1-162 (used in correlation generation)
            double b1,b2,result;            //Temp variables to hold bounds for random number generation when generating correlations

            a = r.Next(1,data.Count());      //Generate random number between 1-162
            if (a != data.Count() - 1)       //If index does not overflow set b1 and b2 using increment
            {
                b1 = data[a];
                b2 = data[a+1];
            }
            else                            //Else set b1 and b2 using decrement
            {
                b1 = data[a-1];
                b2 = data[a];
            }
            result = r.NextDouble() * (b2-b1) + b1;    //Generate random double between b1 and b2
            return result;
        }

        //First Coverage Test Function which tests a 2o3 enviornment with random reliabilities and correlations within the theoretical bounds
        static void CoverageTest1(string system, int runs)
        {
            Random r = new Random();        //Random Number Generator to generate reliabilities
            int success = 0;                //Initialize success counter to zero
            Build(system);                  //Import cuts and paths of system as well as set number of components

            for (int i = 0; i < runs; i++)  //Loop for the appropriate number of test runs 
            {
                Console.WriteLine(i);
                Bm.Clear();                     //Clear B matrix and reinitialize the vector with 2 placeholders
                Bm.Add(new double[] {0});
                Bm.Add(new double[] {0});
                GenerateComponents(0.001, 0.999); //Generate n random reliabilities
                GenerateCorrelations();           //Generate correlations within the theoretical bounds
                GenerateSigma();                  //Generate Sigma values for use in B matrix

                //Clear out list of previous conditional and total probability errors
                ConditionalPErrors.Clear();
                TotalPErrors.Clear();
                //Enumerate Tree
                FindReliability();
                //Check to see if there were any errors when enumerating the tree, if there were no errors increment success
                if (ConditionalPErrors.Count == 0 && TotalPErrors.Count == 0)
                    success++;
            }
            double success_rate = (double)success / (double)runs;
            Console.WriteLine("Coverage Test 1:\n" + success + " successful enumerations after " + runs + " test runs");
            Console.WriteLine("Success Rate: " + success_rate);
        }

        //Second Coverage test which tests a 2o3 enviornment sampling reliabilties from a beta distribution and sampling correlations from previous data set 
        static void CoverageTest2(string system, int runs, bool correction = true)
        {
            bool fail = false;              //Fail flag for non-correction method
            int r,c_index,sum = 0;          //Row index and correlations index used in transposing correaltions into the c array
            double[] correlations;          //Generated correlation values
            double[] c;                     //Temp variable used in correlation matrix construction
            int success = 0;                //Initialize success counter to zero
            Build(system);                  //Import cuts and paths of system as well as set number of components
            //Create Beta Distribution with alpha = 1000 and beta = 2, add 10,000 samples to double array and compute average
            Beta beta = new Beta(1000,2);
            double[] samples = new double[number*runs];
            beta.Samples(samples);
            double average = samples.Average(); 
            //Make sure average of beta distribution is within the bounds
            if (average < 0.997 || average > 0.999)
            {
                Console.WriteLine("Beta Distribution Error: Average sample = " + average.ToString("#.######"));
                return;
            }

            //Correlation data set
            double[] data_set = new double[] {-0.000832051, -0.00077379, -0.000700997, -0.000556345, -0.000551149, -0.000546957, -0.00046025, -0.000377824, -0.000334054, -0.00032533, -0.00030337, -0.000289878, -0.000269457, -0.000250505, -0.000222839, -0.000216432, -0.000192773, -0.000172419, -0.000171405, -0.000150472, -0.00014154, -0.000139644, -0.000136932, -0.00013589, -0.000119422, -0.000118307, -0.000117407, -0.0000959645, -0.0000903688, -0.0000808274, -0.000078077, -0.0000753715, -0.0000717061, -0.0000698334, -0.0000678329, -0.000066352, -0.0000651196, -0.0000573269, -0.0000523427, -0.0000432869, -0.0000420424, -0.0000359503, -0.0000332959, -0.0000292295, -0.0000229813, -0.0000228065, -0.0000168529, -0.0000151666, -0.0000145606, -0.0000139291, -0.0000135653, -0.0000126496, -0.0000111359, -0.00000282844, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0000193941, 0.000506082, 0.00114057, 0.00117882, 0.00171039, 0.00185772, 0.00209051, 0.00232217, 0.00238901, 0.0031334, 0.00338147, 0.00348186, 0.00396074, 0.00484534, 0.00636997, 0.00709894, 0.00750422, 0.011968, 0.0140486, 0.0187122, 0.0384815, 0.0503591, 0.0794512, 0.0862146, 0.0981174, 0.120224, 0.136748, 0.167731, 0.198908, 0.205563, 0.220539, 0.227669, 0.407697, 0.441613, 0.530991, 0.58726};

            for (int i = 0; i < runs; i++)  //Loop for the appropriate number of test runs 
            {
                Console.WriteLine(i);
                //Reinitalize fail flag for non-correction method
                fail = false;
                //Clear B matrix and reinitialize the vector with 2 placeholders
                Bm.Clear();                     
                Bm.Add(new double[] {0});
                Bm.Add(new double[] {0});
                //Clear component reliabilities and correaltions then add n samples as component reliabilities 
                component_reliabilities = new double[number];
                for (int j = 0; j < number; j++)
                    component_reliabilities[j] = samples[number*i+j];
                //Generate correlation bounds
                List<CorrelationPairs> bounds = CorrelationBounds();
                correlations = new double[bounds.Count()];
                //Continue generating correlations using data set as reference until all correlations are within theoretical bounds
                if (correction)
                {
                    for (int j = 0; j < bounds.Count(); j++)
                    {
                        //Generate correaltion value using data set, then repeat if outside theoretical bounds
                        correlations[j] = PseudoSample(data_set);
                        if (correlations[j] > bounds[j].max || correlations[j] < bounds[j].min)
                            j--;
                    }
                }
                //Otherwise generate correlations and fail iteration if a correlation is outside theoretical bounds
                else
                {
                    for (int j = 0; j < bounds.Count(); j++)
                    {
                        correlations[j] = PseudoSample(data_set);
                        if (correlations[j] > bounds[j].max || correlations[j] < bounds[j].min)
                        {
                            fail = true;
                            break;
                        }
                    }
                    if (fail)
                        continue;
                }
                //Enter correlations into a double array then generate matrix
                c = new double[number*number];
                r = 0;
                c_index = 0;
                for (int j = 0; j < number*number;)
                {
                    //Initialize summation
                    sum = r-1;
                    //Enter previous correlations at r > k
                    for (int k = 0; k <  r; k++)
                    {
                        if(k == 0)
                            c[j++] = correlations[r-1];
                        else
                        {
                            c[j++] = correlations[(sum)+(number-k-1)];
                            sum += (number - k - 1);
                        }
                    }
                    //Enter correlation of 1 at r = k
                    c[j++] = 1;
                    //Enter new correlations at r < k
                    while (j % number != 0)
                        c[j++] = correlations[c_index++];
                    //Increment row index
                    r++;
                }
                correlations_matrix = Matrix<double>.Build.Dense(number,number,c);
                GenerateSigma();    //Generate Sigma values for use in B matrix

                //Clear out list of previous conditional and total probability errors
                ConditionalPErrors.Clear();
                TotalPErrors.Clear();
                //Enumerate Tree
                FindReliability();
                //Check to see if there were any errors when enumerating the tree, if there were no errors increment success
                if (ConditionalPErrors.Count == 0 && TotalPErrors.Count == 0)
                    success++;
            }
            double success_rate = (double)success / (double)runs;
            if (correction)
                Console.WriteLine("Coverage Test 2 with Correction:\n" + success + " successful enumerations after " + runs + " test runs");
            else
                Console.WriteLine("Coverage Test 2 without Correction:\n" + success + " successful enumerations after " + runs + " test runs");
            Console.WriteLine("Success Rate: " + success_rate);
        }

        static List<int[]> WolframPerms(List<int[]> perms)
        {
            List<int> temp = new List<int>();
            List<int[]> newList = new List<int[]>();

            foreach(int[] i in perms)
            {
                temp.Clear();
                for(int j = 0; j < i.Length; j++)
                {
                    if(i[j] == 1)
                        temp.Add(j+1);
                }
                newList.Add(temp.ToArray());
            }
            return newList;
        }
        
        //Function which builds the appropriate component network based on input string (does not generate b vectors only intializes)
        static void Build(string s)
        {
            double[] correlations;
            Bm.Clear();
            Bm.Add(new double[] {0});
            Bm.Add(new double[] {0});
            CS.Clear();
            PS.Clear();

            switch (s)
            {
                case "2s":                      //2 component series
                    number = 2;
                    component_reliabilities = new double[] {0.8, 0.75};
                    /* correlations = new double[] { 1.00, -0.28,
                                                -0.28,  1.00}; */
                    correlations = new double[] { 1.0000000, 0.1377770,
                                                  0.1377770, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(2,2,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    PS.Add(new int[] {1,2});
                    break;
                case "2p":                      //2 component parallel
                    number = 2;
                    component_reliabilities = new double[] {0.8, 0.75};
                    /* correlations = new double[] { 1.00, -0.28,
                                                -0.28,  1.00}; */
                    correlations = new double[] { 1.0000000, 0.1377770,
                                                  0.1377770, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(2,2,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    break;
                case "3s":                      //3 component series
                    number = 3;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7};
                    /* correlations = new double[] { 1.00, -0.28, 0.1,
                                                -0.28,  1.00, 0.2,
                                                0.1,   0.2,  1.00}; */
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228,
                                                  0.1377770, 1.0000000, 0.1108590,
                                                  0.0326228, 0.1108590, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(3,3,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    CS.Add(new int[] {3});
                    PS.Add(new int[] {1,2,3});
                    break;
                case "3p":                      //3 component parallel
                    number = 3;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7};
                    /* correlations = new double[] { 1.00, -0.28, 0.1,
                                                -0.28,  1.00, 0.2,
                                                0.1,   0.2,  1.00}; */
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228,
                                                  0.1377770, 1.0000000, 0.1108590,
                                                  0.0326228, 0.1108590, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(3,3,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    PS.Add(new int[] {3});
                    break;
                case "4s":                      //4 component series
                    number = 4;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15,
                                                 -0.28,  1.00,   0.2, 0.30,
                                                   0.1,   0.2,  1.00, 0.32,
                                                 -0.15,  0.30,  0.32, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.0719920,
                                                  0.1377770, 1.0000000, 0.1108590, 0.0981194,
                                                  0.0326228, 0.1108590, 1.0000000, 0.0873844,
                                                  0.0719920, 0.0981194, 0.0873844, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(4,4,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    CS.Add(new int[] {3});
                    CS.Add(new int[] {4});
                    PS.Add(new int[] {1,2,3,4});
                    break;
                case "4p":                      //4 component parallel
                    number = 4;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15,
                                                 -0.28,  1.00,   0.2, 0.30,
                                                   0.1,   0.2,  1.00, 0.32,
                                                 -0.15,  0.30,  0.32, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.0719920,
                                                  0.1377770, 1.0000000, 0.1108590, 0.0981194,
                                                  0.0326228, 0.1108590, 1.0000000, 0.0873844,
                                                  0.0719920, 0.0981194, 0.0873844, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(4,4,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3,4});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    PS.Add(new int[] {3});
                    PS.Add(new int[] {4});
                    break;
                case "5s":                      //5 component series
                    number = 5;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.0719920, -0.0486818,
                                                  0.1377770, 1.0000000, 0.1108590, 0.0981194,  0.1800180,
                                                  0.0326228, 0.1108590, 1.0000000, 0.0873844,  0.1303380,
                                                  0.0719920, 0.0981194, 0.0873844, 1.0000000,  0.0248845,
                                                 -0.0486818, 0.1800180, 0.1303380, 0.0248845,  1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(5,5,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    CS.Add(new int[] {3});
                    CS.Add(new int[] {4});
                    CS.Add(new int[] {5});
                    PS.Add(new int[] {1,2,3,4,5});
                    break;
                case "5p":                      //5 component parallel 
                    number = 5;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.0719920, -0.0486818,
                                                  0.1377770, 1.0000000, 0.1108590, 0.0981194,  0.1800180,
                                                  0.0326228, 0.1108590, 1.0000000, 0.0873844,  0.1303380,
                                                  0.0719920, 0.0981194, 0.0873844, 1.0000000,  0.0248845,
                                                 -0.0486818, 0.1800180, 0.1303380, 0.0248845,  1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(5,5,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3,4,5});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    PS.Add(new int[] {3});
                    PS.Add(new int[] {4});
                    PS.Add(new int[] {5});
                    break;
                case "6s":                      //6 component series 
                    number = 6;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, 
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, 
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, 
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, 
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000};
                    correlations_matrix = Matrix<double>.Build.Dense(6,6,correlations);
                    GenerateSigma();
                    PS.Add(new int[] {1,2,3,4,5,6});
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    CS.Add(new int[] {3});
                    CS.Add(new int[] {4});
                    CS.Add(new int[] {5});
                    CS.Add(new int[] {6});
                    break;
                case "6p":                      //6 component parallel 
                    number = 6;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, 
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, 
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, 
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, 
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000};
                    correlations_matrix = Matrix<double>.Build.Dense(6,6,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3,4,5,6});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    PS.Add(new int[] {3});
                    PS.Add(new int[] {4});
                    PS.Add(new int[] {5});
                    PS.Add(new int[] {6});
                    break;
                case "7s":                      //7 component series 
                    number = 7;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77, 0.67};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, -0.066384200, 
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, -0.059739400, 
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, -0.056192500, 
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, -0.000253982, 
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747,  0.127017000, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000,  0.014421100, 
                                                 -0.0663842,-0.0597394,-0.0561925,-0.000253982,  0.12701700, 0.01442110,  1.000000000};
                    correlations_matrix = Matrix<double>.Build.Dense(7,7,correlations);
                    GenerateSigma();
                    PS.Add(new int[] {1,2,3,4,5,6,7});
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    CS.Add(new int[] {3});
                    CS.Add(new int[] {4});
                    CS.Add(new int[] {5});
                    CS.Add(new int[] {6});
                    CS.Add(new int[] {7});
                    break;
                case "7p":                      //7 component parallel 
                    number = 7;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77, 0.67};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, -0.066384200, 
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, -0.059739400, 
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, -0.056192500, 
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, -0.000253982, 
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747,  0.127017000, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000,  0.014421100, 
                                                 -0.0663842,-0.0597394,-0.0561925,-0.000253982,  0.12701700, 0.01442110,  1.000000000};
                    correlations_matrix = Matrix<double>.Build.Dense(7,7,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3,4,5,6,7});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    PS.Add(new int[] {3});
                    PS.Add(new int[] {4});
                    PS.Add(new int[] {5});
                    PS.Add(new int[] {6});
                    PS.Add(new int[] {7});
                    break;
                case "8s":                      //8 component series 
                    number = 8;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77, 0.67,0.5};
                    /* correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23, 0.22,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,-0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17, 0.07,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 0.25,
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16, 0.10,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,-0.07,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00, 0.35,
                                                  0.22, -0.16,  0.07, 0.25, 0.10,-0.07, 0.35, 1.00}; */
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, -0.066384200, 0.0178629,
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, -0.059739400,-0.0346766,
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, -0.056192500,-0.0402154,
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, -0.000253982, 0.0286271,
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747,  0.127017000,-0.0166433, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000,  0.014421100,-0.0140690, 
                                                 -0.0663842,-0.0597394,-0.0561925,-0.000253982,  0.12701700, 0.01442110,  1.000000000,-0.0305813,
                                                  0.0178629,-0.0346766,-0.0402154, 0.028627100, -0.01664330, -0.0140690, -0.030581300, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(8,8,correlations);
                    GenerateSigma();
                    PS.Add(new int[] {1,2,3,4,5,6,7,8});
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    CS.Add(new int[] {3});
                    CS.Add(new int[] {4});
                    CS.Add(new int[] {5});
                    CS.Add(new int[] {6});
                    CS.Add(new int[] {7});
                    CS.Add(new int[] {8});
                    break;
                case "8p":                      //8 component parallel 
                    number = 8;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77, 0.67,0.5};
                    /* correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23, 0.22,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,-0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17, 0.07,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 0.25,
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16, 0.10,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,-0.07,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00, 0.35,
                                                  0.22, -0.16,  0.07, 0.25, 0.10,-0.07, 0.35, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, -0.066384200, 0.0178629,
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, -0.059739400,-0.0346766,
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, -0.056192500,-0.0402154,
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, -0.000253982, 0.0286271,
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747,  0.127017000,-0.0166433, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000,  0.014421100,-0.0140690, 
                                                 -0.0663842,-0.0597394,-0.0561925,-0.000253982,  0.12701700, 0.01442110,  1.000000000,-0.0305813,
                                                  0.0178629,-0.0346766,-0.0402154, 0.028627100, -0.01664330, -0.0140690, -0.030581300, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(8,8,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3,4,5,6,7,8});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    PS.Add(new int[] {3});
                    PS.Add(new int[] {4});
                    PS.Add(new int[] {5});
                    PS.Add(new int[] {6});
                    PS.Add(new int[] {7});
                    PS.Add(new int[] {8});
                    break;
                case "5bn":                     //5 component bridge network
                    number = 5;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.0719920, -0.0486818,
                                                  0.1377770, 1.0000000, 0.1108590, 0.0981194,  0.1800180,
                                                  0.0326228, 0.1108590, 1.0000000, 0.0873844,  0.1303380,
                                                  0.0719920, 0.0981194, 0.0873844, 1.0000000,  0.0248845,
                                                 -0.0486818, 0.1800180, 0.1303380, 0.0248845,  1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(5,5,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2});
                    CS.Add(new int[] {4,5});
                    CS.Add(new int[] {1,3,5});
                    CS.Add(new int[] {2,3,4});
                    PS.Add(new int[] {1,4});
                    PS.Add(new int[] {2,3,4});
                    PS.Add(new int[] {2,5});
                    PS.Add(new int[] {1,3,5});
                    break;
                case "10bn":                     //10 component bridge network
                    number = 10;
                    component_reliabilities = new double[] {0.97499, 0.99267, 0.97783, 0.95330, 0.95391, 0.95755, 0.97215, 0.98216, 0.99780, 0.99689};
                    correlations = new double[] {1.00000, 0.03505, 0.06065, 0.04038, 0.06685, 0.03315, 0.06523, 0.04131, -0.00012, 0.01537,
                                                    0.03505, 1.00000, 0.04604, 0.00359, 0.03011, 0.02811, 0.04971, 0.01770, 0.05030,
                                                    0.01026, 0.06065, 0.04604, 1.00000, 0.01324, 0.04172, 0.05361, -0.00036, 0.04143,
                                                    0.03017, 0.02362, 0.04038, 0.00359, 0.01324, 1.00000, 0.02539, 0.05066, 0.01900,
                                                    0.05321, 0.00434, 0.01259, 0.06685, 0.03011, 0.04172, 0.02539, 1.00000, 0.09488,
                                                    0.01285, 0.05671, 0.00744, -0.00002, 0.03315, 0.02811, 0.05361, 0.05066, 0.09488,
                                                    1.00000, 0.00497, 0.00014, 0.01566, 0.01765, 0.06523, 0.04971, -0.00036, 0.01900,
                                                    0.01285, 0.00497, 1.00000, 0.02312, 0.01378, 0.01076, 0.04131, 0.01770, 0.04143,
                                                    0.05321, 0.05671, 0.00014, 0.02312, 1.00000, 0.02180, 0.02158, -0.00012, 0.05030,
                                                    0.03017, 0.00434, 0.00744, 0.01566, 0.01378, 0.02180, 1.00000, 0.08048, 0.01537,
                                                    0.01026, 0.02362, 0.01259, -0.00002, 0.01765, 0.01076, 0.02158, 0.08048, 1.00000};
                    correlations_matrix = Matrix<double>.Build.Dense(10,10,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,3}); CS.Add(new int[] {1,4}); CS.Add(new int[] {2,3}); CS.Add(new int[] {2,4});
                    CS.Add(new int[] {7,9}); CS.Add(new int[] {7,10}); CS.Add(new int[] {8,9}); CS.Add(new int[] {8,10});
                    CS.Add(new int[] {1,5,9}); CS.Add(new int[] {1,6,9}); CS.Add(new int[] {1,5,10}); CS.Add(new int[] {1,6,10});
                    CS.Add(new int[] {2,5,9}); CS.Add(new int[] {2,6,9}); CS.Add(new int[] {2,5,10}); CS.Add(new int[] {2,6,10});
                    CS.Add(new int[] {3,6,7}); CS.Add(new int[] {3,5,7}); CS.Add(new int[] {3,6,8}); CS.Add(new int[] {3,5,8});
                    CS.Add(new int[] {4,6,7}); CS.Add(new int[] {4,5,7}); CS.Add(new int[] {4,6,8}); CS.Add(new int[] {4,5,8});

                    PS.Add(new int[] {1,2,7,8});
                    PS.Add(new int[] {3,4,9,10});
                    PS.Add(new int[] {1,2,5,6,9,10});
                    PS.Add(new int[] {3,4,6,5,7,8});
                    break;
                case "20bn":                     //20 component bridge network
                    number = 20;
                    component_reliabilities = new double[] {0.96166, 0.98371, 0.96786, 0.97608, 0.95177, 0.96956, 0.99830, 0.97229, 
                    0.96720, 0.98829, 0.97603, 0.96054, 0.98679, 0.99279, 0.97152, 0.95646, 0.96767, 0.98632, 0.97210, 0.99855};
                    correlations = new double[] {1.00000, 0.01477, 0.04543, -0.00060, 0.01613, 0.04260, 0.00018, 0.01548, 0.02258, 0.00666, 0.00835, 0.01386, 0.00335, 0.01451, 0.04165, 0.04061, 0.01109, 0.02162, 0.01295, 0.00284, 0.01477,
                                                    1.00000, 0.03031, 0.00887, 0.01092, -0.00113, 0.01444, 0.01415, 0.00724, 0.02365, 0.00455, 0.02643, 0.04128, 0.01151, 0.02741, 0.00375, 0.00777, 0.02480, 0.03120, 0.00376, 0.04543,
                                                    0.03031, 1.00000, 0.02775, 0.01510, 0.02799, 0.00847, 0.03712, 0.01322, 0.00403, 0.01615, 0.03089, 0.02642, 0.01666, 0.01649, 0.02856, 0.02703, 0.02585, 0.04178, 0.00318, -0.00060,
                                                    0.00887, 0.02775, 1.00000, 0.00549, 0.00168, 0.01311, 0.00445, -0.00119, 0.02509, 0.02376, 0.00881, 0.02446, 0.00919, 0.04029, 0.02029, 0.04226, 0.01627, 0.01878, 0.00029, 0.01613,
                                                    0.01092, 0.01510, 0.00549, 1.00000, -0.00179, 0.00036, 0.03469, 0.02953, 0.00148, 0.01384, 0.00341, 0.02131, 0.00912, 0.02561, 0.01730, 0.03398, 0.01956, 0.00426, 0.00279, 0.04260,
                                                    -0.00113, 0.02799, 0.00168, -0.00179, 1.00000, 0.00626, 0.03953, 0.03992, 0.00054, 0.02060, 0.03603, 0.00387, 0.01541, 0.03211, 0.00250, 0.02332, 0.00368, 0.04473, 0.00274, 0.00018,
                                                    0.01444, 0.00847, 0.01311, 0.00036, 0.00626, 1.00000, 0.00999, 0.00787, 0.00728, 0.00370, 0.00302, 0.01201, 0.01131, 0.01053, 0.00124, 0.00219, 0.01346, 0.01151, 0.00722, 0.01548,
                                                    0.01415, 0.03712, 0.00445, 0.03469, 0.03953, 0.00999, 1.00000, 0.02754, 0.00958, 0.01456, 0.02618, 0.01659, 0.01263, 0.04884, 0.01196, 0.03365, 0.00445, 0.01899, 0.01059, 0.02258,
                                                    0.00724, 0.01322, -0.00119, 0.02953, 0.03992, 0.00787, 0.02754, 1.00000, 0.00560, 0.03636, 0.01867, 0.00718, 0.01439, 0.02768, 0.01107, 0.00769, 0.01741, 0.03096, 0.01009, 0.00666,
                                                    0.02365, 0.00403, 0.02509, 0.00148, 0.00054, 0.00728, 0.00958, 0.00560, 1.00000, -0.00071, 0.01417, 0.02012, 0.01846, 0.02824, 0.01682, 0.01957, 0.02819, 0.00820, 0.00927, 0.00835,
                                                    0.00455, 0.01615, 0.02376, 0.01384, 0.02060, 0.00370, 0.01456, 0.03636, -0.00071, 1.00000, 0.02551, 0.01919, 0.00867, 0.01478, 0.01590, 0.04284, 0.02309, 0.04371, 0.01147, 0.01386,
                                                    0.02643, 0.03089, 0.00881, 0.00341, 0.03603, 0.00302, 0.02618, 0.01867, 0.01417, 0.02551, 1.00000, 0.00137, 0.01141, 0.04193, 0.01561, 0.00326, 0.02356, -0.00011, 0.00571, 0.00335,
                                                    0.04128, 0.02642, 0.02446, 0.02131, 0.00387, 0.01201, 0.01659, 0.00718, 0.02012, 0.01919, 0.00137, 1.00000, 0.01938, 0.02028, 0.01245, 0.02394, 0.02911, 0.03169, 0.00545, 0.01451,
                                                    0.01151, 0.01666, 0.00919, 0.00912, 0.01541, 0.01131, 0.01263, 0.01439, 0.01846, 0.00867, 0.01141, 0.01938, 1.00000, 0.00769, 0.00112, 0.00074, 0.00083, 0.00185, 0.01466, 0.04165,
                                                    0.02741, 0.01649, 0.04029, 0.02561, 0.03211, 0.01053, 0.04884, 0.02768, 0.02824, 0.01478, 0.04193, 0.02028, 0.00769, 1.00000, 0.00107, 0.02440, 0.00333, 0.01920, 0.00588, 0.04061,
                                                    0.00375, 0.02856, 0.02029, 0.01730, 0.00250, 0.00124, 0.01196, 0.01107, 0.01682, 0.01590, 0.01561, 0.01245, 0.00112, 0.00107, 1.00000, 0.01275, 0.00160, 0.03603, 0.00124, 0.01109,
                                                    0.00777, 0.02703, 0.04226, 0.03398, 0.02332, 0.00219, 0.03365, 0.00769, 0.01957, 0.04284, 0.00326, 0.02394, 0.00074, 0.02440, 0.01275, 1.00000, 0.01921, 0.03032, 0.00219, 0.02162,
                                                    0.02480, 0.02585, 0.01627, 0.01956, 0.00368, 0.01346, 0.00445, 0.01741, 0.02819, 0.02309, 0.02356, 0.02911, 0.00083, 0.00333, 0.00160, 0.01921, 1.00000, 0.02096, 0.00958, 0.01295,
                                                    0.03120, 0.04178, 0.01878, 0.00426, 0.04473, 0.01151, 0.01899, 0.03096, 0.00820, 0.04371, -0.00011, 0.03169, 0.00185, 0.01920, 0.03603, 0.03032, 0.02096, 1.00000, 0.00449, 0.00284,
                                                    0.00376, 0.00318, 0.00029, 0.00279, 0.00274, 0.00722, 0.01059, 0.01009, 0.00927, 0.01147, 0.00571, 0.00545, 0.01466, 0.00588, 0.00124, 0.00219, 0.00958, 0.00449, 1.00000};
                    correlations_matrix = Matrix<double>.Build.Dense(20,20,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,5}); CS.Add(new int[] {1,6}); CS.Add(new int[] {1,7}); CS.Add(new int[] {1,8});
                    CS.Add(new int[] {2,5}); CS.Add(new int[] {2,6}); CS.Add(new int[] {2,7}); CS.Add(new int[] {2,8});
                    CS.Add(new int[] {3,5}); CS.Add(new int[] {3,6}); CS.Add(new int[] {3,7}); CS.Add(new int[] {3,8});
                    CS.Add(new int[] {4,5}); CS.Add(new int[] {4,6}); CS.Add(new int[] {4,7}); CS.Add(new int[] {4,8});
                    CS.Add(new int[] {16,17}); CS.Add(new int[] {16,18}); CS.Add(new int[] {16,19}); CS.Add(new int[] {16,20});
                    CS.Add(new int[] {15,17}); CS.Add(new int[] {15,18}); CS.Add(new int[] {15,19}); CS.Add(new int[] {15,20});
                    CS.Add(new int[] {14,17}); CS.Add(new int[] {14,18}); CS.Add(new int[] {14,19}); CS.Add(new int[] {14,20});
                    CS.Add(new int[] {13,17}); CS.Add(new int[] {13,18}); CS.Add(new int[] {13,19}); CS.Add(new int[] {13,20});
                    CS.Add(new int[] {1,9,17}); CS.Add(new int[] {1,9,18}); CS.Add(new int[] {1,9,19}); CS.Add(new int[] {1,9,20});
                    CS.Add(new int[] {1,10,17}); CS.Add(new int[] {1,10,18}); CS.Add(new int[] {1,10,19}); CS.Add(new int[] {1,10,20});
                    CS.Add(new int[] {1,11,17}); CS.Add(new int[] {1,11,18}); CS.Add(new int[] {1,11,19}); CS.Add(new int[] {1,11,20});
                    CS.Add(new int[] {1,12,17}); CS.Add(new int[] {1,12,18}); CS.Add(new int[] {1,12,19}); CS.Add(new int[] {1,12,20});
                    CS.Add(new int[] {5,9,13}); CS.Add(new int[] {5,9,14}); CS.Add(new int[] {5,9,15}); CS.Add(new int[] {5,9,16});
                    CS.Add(new int[] {5,10,13}); CS.Add(new int[] {5,10,14}); CS.Add(new int[] {5,10,15}); CS.Add(new int[] {5,10,16});
                    CS.Add(new int[] {5,11,13}); CS.Add(new int[] {5,11,14}); CS.Add(new int[] {5,11,15}); CS.Add(new int[] {5,11,16});
                    CS.Add(new int[] {5,12,13}); CS.Add(new int[] {5,12,14}); CS.Add(new int[] {5,12,15}); CS.Add(new int[] {5,12,16});
                    CS.Add(new int[] {2,9,17}); CS.Add(new int[] {2,9,18}); CS.Add(new int[] {2,9,19}); CS.Add(new int[] {2,9,20});
                    CS.Add(new int[] {2,10,17}); CS.Add(new int[] {2,10,18}); CS.Add(new int[] {2,10,19}); CS.Add(new int[] {2,10,20});
                    CS.Add(new int[] {2,11,17}); CS.Add(new int[] {2,11,18}); CS.Add(new int[] {2,11,19}); CS.Add(new int[] {2,11,20});
                    CS.Add(new int[] {2,12,17}); CS.Add(new int[] {2,12,18}); CS.Add(new int[] {2,12,19}); CS.Add(new int[] {2,12,20});
                    CS.Add(new int[] {6,9,13}); CS.Add(new int[] {6,9,14}); CS.Add(new int[] {6,9,15}); CS.Add(new int[] {6,9,16});
                    CS.Add(new int[] {6,10,13}); CS.Add(new int[] {6,10,14}); CS.Add(new int[] {6,10,15}); CS.Add(new int[] {6,10,16});
                    CS.Add(new int[] {6,11,13}); CS.Add(new int[] {6,11,14}); CS.Add(new int[] {6,11,15}); CS.Add(new int[] {6,11,16});
                    CS.Add(new int[] {6,12,13}); CS.Add(new int[] {6,12,14}); CS.Add(new int[] {6,12,15}); CS.Add(new int[] {6,12,16});
                    CS.Add(new int[] {3,9,17}); CS.Add(new int[] {3,9,18}); CS.Add(new int[] {3,9,19}); CS.Add(new int[] {3,9,20});
                    CS.Add(new int[] {3,10,17}); CS.Add(new int[] {3,10,18}); CS.Add(new int[] {3,10,19}); CS.Add(new int[] {3,10,20});
                    CS.Add(new int[] {3,11,17}); CS.Add(new int[] {3,11,18}); CS.Add(new int[] {3,11,19}); CS.Add(new int[] {3,11,20});
                    CS.Add(new int[] {3,12,17}); CS.Add(new int[] {3,12,18}); CS.Add(new int[] {3,12,19}); CS.Add(new int[] {3,12,20});
                    CS.Add(new int[] {4,9,17}); CS.Add(new int[] {4,9,18}); CS.Add(new int[] {4,9,19}); CS.Add(new int[] {4,9,20});
                    CS.Add(new int[] {4,10,17}); CS.Add(new int[] {4,10,18}); CS.Add(new int[] {4,10,19}); CS.Add(new int[] {4,10,20});
                    CS.Add(new int[] {4,11,17}); CS.Add(new int[] {4,11,18}); CS.Add(new int[] {4,11,19}); CS.Add(new int[] {4,11,20});
                    CS.Add(new int[] {4,12,17}); CS.Add(new int[] {4,12,18}); CS.Add(new int[] {4,12,19}); CS.Add(new int[] {4,12,20});
                    CS.Add(new int[] {7,9,13}); CS.Add(new int[] {7,9,14}); CS.Add(new int[] {7,9,15}); CS.Add(new int[] {7,9,16});
                    CS.Add(new int[] {7,10,13}); CS.Add(new int[] {7,10,14}); CS.Add(new int[] {7,10,15}); CS.Add(new int[] {7,10,16});
                    CS.Add(new int[] {7,11,13}); CS.Add(new int[] {7,11,14}); CS.Add(new int[] {7,11,15}); CS.Add(new int[] {7,11,16});
                    CS.Add(new int[] {7,12,13}); CS.Add(new int[] {7,12,14}); CS.Add(new int[] {7,12,15}); CS.Add(new int[] {7,12,16});
                    CS.Add(new int[] {8,9,13}); CS.Add(new int[] {8,9,14}); CS.Add(new int[] {8,9,15}); CS.Add(new int[] {8,9,16});
                    CS.Add(new int[] {8,10,13}); CS.Add(new int[] {8,10,14}); CS.Add(new int[] {8,10,15}); CS.Add(new int[] {8,10,16});
                    CS.Add(new int[] {8,11,13}); CS.Add(new int[] {8,11,14}); CS.Add(new int[] {8,11,15}); CS.Add(new int[] {8,11,16});
                    CS.Add(new int[] {8,12,13}); CS.Add(new int[] {8,12,14}); CS.Add(new int[] {8,12,15}); CS.Add(new int[] {8,12,16});
                    
                    PS.Add(new int[] {1,2,3,4,13,14,15,16});
                    PS.Add(new int[] {5,6,7,8,17,18,19,20});
                    PS.Add(new int[] {1,2,3,4,9,10,11,12,17,18,19,20});
                    PS.Add(new int[] {5,6,7,8,12,11,10,9,13,14,15,16});
                    break;
                case "22bn":
                    number = 22;
                    component_reliabilities = new double[] {0.98671, 0.88780, 0.86846, 0.99231, 0.97476, 0.89183, 0.97980, 0.96316, 0.92800, 0.96436, 0.94459, 0.87680, 0.91635, 0.91582, 0.87922, 0.94138, 0.98290, 0.88984, 0.89236, 0.92153, 0.91587, 0.96817};
                    correlations = new double [] {1.00000, -0.00025, 0.00105, 0.00169, 0.00447, 0.00174, 0.00568, 0.00415, 0.00210, 0.00452, 0.00172, 0.00302, 0.00064, 0.00354, -0.00029, 0.00453, 0.00511, 0.00077, 0.00143, 0.00368, 0.00116, 0.00107, -0.00025,
                                                    1.00000, 0.00016, 0.00023, 0.00430, 0.00403, 0.00366, 0.00448, 0.00133, 0.00087, 0.00397, 0.00501, -0.00100, 0.00457, 0.00329, 0.00558, 0.00030, 0.00131, 0.00605, 0.00323, 0.00614, 0.00176, 0.00105,
                                                    0.00016, 1.00000, 0.00214, -0.00048, 0.00015, 0.00355, 0.00218, -0.00042, 0.00140, 0.00408, 0.00210, 0.00452, -0.00019, 0.00801, 0.00544, 0.00001, 0.00798, 0.00113, 0.00126, 0.00019, 0.00282, 0.00169,
                                                    0.00023, 0.00214, 1.00000, 0.00476, 0.00096, 0.00165, 0.00043, 0.00250, 0.00444, 0.00228, 0.00103, 0.00044, 0.00225, 0.00064, 0.00223, 0.00007, 0.00020, 0.00237, 0.00052, 0.00010, 0.00009, 0.00447,
                                                    0.00430, -0.00048, 0.00476, 1.00000, 0.00455, 0.00864, 0.00684, 0.00240, 0.00457, 0.00618, 0.00425, 0.00228, 0.00035, 0.00296, 0.00156, 0.00752, 0.00427, 0.00438, 0.00258, 0.00418, 0.00774, 0.00174,
                                                    0.00403, 0.00015, 0.00096, 0.00455, 1.00000, 0.00399, 0.00354, 0.00713, 0.00143, 0.00693, 0.00782, 0.00452, 0.00047, 0.00841, 0.00686, 0.00371, -0.00030, 0.00642, 0.00768, 0.00770, 0.00394, 0.00568,
                                                    0.00366, 0.00355, 0.00165, 0.00864, 0.00399, 1.00000, 0.00118, 0.00122, 0.00115, 0.00480, 0.00128, 0.00416, 0.00138, 0.00007, 0.00166, 0.00756, 0.00375, 0.00052, 0.00119, 0.00397, 0.00307, 0.00415,
                                                    0.00448, 0.00218, 0.00043, 0.00684, 0.00354, 0.00118, 1.00000, 0.00055, 0.00189, 0.00212, 0.00505, 0.00342, 0.00090, 0.00354, 0.00775, 0.00050, 0.00096, 0.00122, 0.00419, 0.00625, -0.00026, 0.00210,
                                                    0.00133, -0.00042, 0.00250, 0.00240, 0.00713, 0.00122, 0.00055, 1.00000, 0.00239, 0.00366, 0.00569, 0.00308, 0.00053, 0.00538, 0.00569, 0.00130, 0.00779, 0.00316, 0.00007, 0.00236, 0.00487, 0.00452,
                                                    0.00087, 0.00140, 0.00444, 0.00457, 0.00143, 0.00115, 0.00189, 0.00239, 1.00000, -0.00040, 0.00028, 0.00590, 0.00185, 0.00433, -0.00002, 0.00370, 0.00447, 0.00467, 0.00215, 0.00602, 0.00288, 0.00172,
                                                    0.00397, 0.00408, 0.00228, 0.00618, 0.00693, 0.00480, 0.00212, 0.00366, -0.00040, 1.00000, 0.00330, 0.00347, 0.00189, -0.00033, 0.00512, 0.00035, 0.00093, 0.00213, 0.00044, 0.00745, 0.00244, 0.00302,
                                                    0.00501, 0.00210, 0.00103, 0.00425, 0.00782, 0.00128, 0.00505, 0.00569, 0.00028, 0.00330, 1.00000, 0.00009, 0.00653, 0.00046, -0.00017, 0.00250, 0.00745, 0.00120, 0.00639, 0.00704, 0.00069, 0.00064,
                                                    -0.00100, 0.00452, 0.00044, 0.00228, 0.00452, 0.00416, 0.00342, 0.00308, 0.00590, 0.00347, 0.00009, 1.00000, 0.00393, 0.00539, 0.00663, 0.00043, 0.00774, 0.00330, 0.00542, 0.00361, 0.00058, 0.00354,
                                                    0.00457, -0.00019, 0.00225, 0.00035, 0.00047, 0.00138, 0.00090, 0.00053, 0.00185, 0.00189, 0.00653, 0.00393, 1.00000, 0.00383, -0.00001, 0.00360, -0.00048, 0.00449, 0.00195, 0.00301, 0.00108, -0.00029,
                                                    0.00329, 0.00801, 0.00064, 0.00296, 0.00841, 0.00007, 0.00354, 0.00538, 0.00433, -0.00033, 0.00046, 0.00539, 0.00383, 1.00000, 0.00237, 0.00090, 0.00792, 0.00872, 0.00465, -0.00047, 0.00382, 0.00453,
                                                    0.00558, 0.00544, 0.00223, 0.00156, 0.00686, 0.00166, 0.00775, 0.00569, -0.00002, 0.00512, -0.00017, 0.00663, -0.00001, 0.00237, 1.00000, 0.00416, 0.00327, 0.00470, 0.00814, 0.00652, -0.00031, 0.00511,
                                                    0.00030, 0.00001, 0.00007, 0.00752, 0.00371, 0.00756, 0.00050, 0.00130, 0.00370, 0.00035, 0.00250, 0.00043, 0.00360, 0.00090, 0.00416, 1.00000, 0.00217, 0.00182, -0.00006, -0.00001, 0.00070, 0.00077,
                                                    0.00131, 0.00798, 0.00020, 0.00427, -0.00030, 0.00375, 0.00096, 0.00779, 0.00447, 0.00093, 0.00745, 0.00774, -0.00048, 0.00792, 0.00327, 0.00217, 1.00000, 0.00741, 0.00346, 0.00766, 0.00308, 0.00143,
                                                    0.00605, 0.00113, 0.00237, 0.00438, 0.00642, 0.00052, 0.00122, 0.00316, 0.00467, 0.00213, 0.00120, 0.00330, 0.00449, 0.00872, 0.00470, 0.00182, 0.00741, 1.00000, 0.00529, 0.00513, 0.00501, 0.00368,
                                                    0.00323, 0.00126, 0.00052, 0.00258, 0.00768, 0.00119, 0.00419, 0.00007, 0.00215, 0.00044, 0.00639, 0.00542, 0.00195, 0.00465, 0.00814, -0.00006, 0.00346, 0.00529, 1.00000, 0.00916, 0.00220, 0.00116,
                                                    0.00614, 0.00019, 0.00010, 0.00418, 0.00770, 0.00397, 0.00625, 0.00236, 0.00602, 0.00745, 0.00704, 0.00361, 0.00301, -0.00047, 0.00652, -0.00001, 0.00766, 0.00513, 0.00916, 1.00000, 0.00129, 0.00107,
                                                    0.00176, 0.00282, 0.00009, 0.00774, 0.00394, 0.00307, -0.00026, 0.00487, 0.00288, 0.00244, 0.00069, 0.00058, 0.00108, 0.00382, -0.00031, 0.00070, 0.00308, 0.00501, 0.00220, 0.00129, 1.00000};
                    correlations_matrix = Matrix<double>.Build.Dense(22,22,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1, 16});
                    CS.Add(new int[] {1, 17}); CS.Add(new int[] {2, 16}); CS.Add(new int[] {2, 17}); CS.Add(new int[] {7, 22});                  CS.Add(new int[] {7, 21});
                    CS.Add(new int[] {7, 20}); CS.Add(new int[] {6, 22}); CS.Add(new int[] {6, 21}); CS.Add(new int[] {6, 20});
                    CS.Add(new int[] {5, 22});
                    CS.Add(new int[] {5, 21});
                    CS.Add(new int[] {5, 20});

                    CS.Add(new int[] {1, 9, 18});
                    CS.Add(new int[] {1, 8, 18});

                    CS.Add(new int[] {2, 9, 18});
                    CS.Add(new int[] {2, 8, 18});

                    CS.Add(new int[] {1, 9, 19});
                    CS.Add(new int[] {1, 8, 19});

                    CS.Add(new int[] {2, 9, 19});
                    CS.Add(new int[] {2, 8, 19});

                    CS.Add(new int[] {1, 8, 13,15,20});
                    CS.Add(new int[] {1, 8, 12,15,20});
                    CS.Add(new int[] {1, 8, 11,15,20});
                    CS.Add(new int[] {1, 8, 10,15,20});

                    CS.Add(new int[] {1, 8, 13,14,20});
                    CS.Add(new int[] {1, 8, 12,14,20});
                    CS.Add(new int[] {1, 8, 11,14,20});
                    CS.Add(new int[] {1, 8, 10,14,20});

                    CS.Add(new int[] {1, 8, 13,15,21});
                    CS.Add(new int[] {1, 8, 12,15,21});
                    CS.Add(new int[] {1, 8, 11,15,21});
                    CS.Add(new int[] {1, 8, 10,15,21});

                    CS.Add(new int[] {1, 8, 13,14,21});
                    CS.Add(new int[] {1, 8, 12,14,21});
                    CS.Add(new int[] {1, 8, 11,14,21});
                    CS.Add(new int[] {1, 8, 10,14,21});

                    CS.Add(new int[] {1, 8, 13,15,22});
                    CS.Add(new int[] {1, 8, 12,15,22});
                    CS.Add(new int[] {1, 8, 11,15,22});
                    CS.Add(new int[] {1, 8, 10,15,22});

                    CS.Add(new int[] {1, 8, 13,14,22});
                    CS.Add(new int[] {1, 8, 12,14,22});
                    CS.Add(new int[] {1, 8, 11,14,22});
                    CS.Add(new int[] {1, 8, 10,14,22});

                    CS.Add(new int[] {2, 8, 13,15,20});
                    CS.Add(new int[] {2, 8, 12,15,20});
                    CS.Add(new int[] {2, 8, 11,15,20});
                    CS.Add(new int[] {2, 8, 10,15,20});

                    CS.Add(new int[] {2, 8, 13,14,20});
                    CS.Add(new int[] {2, 8, 12,14,20});
                    CS.Add(new int[] {2, 8, 11,14,20});
                    CS.Add(new int[] {2, 8, 10,14,20});

                    CS.Add(new int[] {2, 8, 13,15,21});
                    CS.Add(new int[] {2, 8, 12,15,21});
                    CS.Add(new int[] {2, 8, 11,15,21});
                    CS.Add(new int[] {2, 8, 10,15,21});

                    CS.Add(new int[] {2, 8, 13,14,21});
                    CS.Add(new int[] {2, 8, 12,14,21});
                    CS.Add(new int[] {2, 8, 11,14,21});
                    CS.Add(new int[] {2, 8, 10,14,21});

                    CS.Add(new int[] {2, 8, 13,15,22});
                    CS.Add(new int[] {2, 8, 12,15,22});
                    CS.Add(new int[] {2, 8, 11,15,22});
                    CS.Add(new int[] {2, 8, 10,15,22});

                    CS.Add(new int[] {2, 8, 13,14,22});
                    CS.Add(new int[] {2, 8, 12,14,22});
                    CS.Add(new int[] {2, 8, 11,14,22});
                    CS.Add(new int[] {2, 8, 10,14,22});

                    CS.Add(new int[] {2, 9, 13,15,20});
                    CS.Add(new int[] {2, 9, 12,15,20});
                    CS.Add(new int[] {2, 9, 11,15,20});
                    CS.Add(new int[] {2, 9, 10,15,20});

                    CS.Add(new int[] {2, 9, 13,14,20});
                    CS.Add(new int[] {2, 9, 12,14,20});
                    CS.Add(new int[] {2, 9, 11,14,20});
                    CS.Add(new int[] {2, 9, 10,14,20});

                    CS.Add(new int[] {2, 9, 13,15,21});
                    CS.Add(new int[] {2, 9, 12,15,21});
                    CS.Add(new int[] {2, 9, 11,15,21});
                    CS.Add(new int[] {2, 9, 10,15,21});

                    CS.Add(new int[] {2, 9, 13,14,21});
                    CS.Add(new int[] {2, 9, 12,14,21});
                    CS.Add(new int[] {2, 9, 11,14,21});
                    CS.Add(new int[] {2, 9, 10,14,21});

                    CS.Add(new int[] {2, 9, 13,15,22});
                    CS.Add(new int[] {2, 9, 12,15,22});
                    CS.Add(new int[] {2, 9, 11,15,22});
                    CS.Add(new int[] {2, 9, 10,15,22});

                    CS.Add(new int[] {2, 9, 13,14,22});
                    CS.Add(new int[] {2,9, 12,14,22});
                    CS.Add(new int[] {2, 9, 11,14,22});
                    CS.Add(new int[] {2, 9, 10,14,22});

                    CS.Add(new int[] {1, 9, 13,15,20});
                    CS.Add(new int[] {1, 9, 12,15,20});
                    CS.Add(new int[] {1, 9, 11,15,20});
                    CS.Add(new int[] {1, 9, 10,15,20});

                    CS.Add(new int[] {1, 9, 13,14,20});
                    CS.Add(new int[] {1, 9, 12,14,20});
                    CS.Add(new int[] {1, 9, 11,14,20});
                    CS.Add(new int[] {1, 9, 10,14,20});

                    CS.Add(new int[] {1, 9, 13,15,21});
                    CS.Add(new int[] {1, 9, 12,15,21});
                    CS.Add(new int[] {1, 9, 11,15,21});
                    CS.Add(new int[] {1, 9, 10,15,21});

                    CS.Add(new int[] {1, 9, 13,14,21});
                    CS.Add(new int[] {1, 9, 12,14,21});
                    CS.Add(new int[] {1, 9, 11,14,21});
                    CS.Add(new int[] {1, 9, 10,14,21});

                    CS.Add(new int[] {1, 9, 13,15,22});
                    CS.Add(new int[] {1, 9, 12,15,22});
                    CS.Add(new int[] {1, 9, 11,15,22});
                    CS.Add(new int[] {1, 9, 10,15,22});

                    CS.Add(new int[] {1, 9, 13,14,22});
                    CS.Add(new int[] {1,9, 12,14,22});
                    CS.Add(new int[] {1, 9, 11,14,22});
                    CS.Add(new int[] {1, 9, 10,14,22});
                    CS.Add(new int[] {3,14, 22});
                    CS.Add(new int[] {3,14, 21});
                    CS.Add(new int[] {3,14, 20});

                    CS.Add(new int[] {4,14, 22});
                    CS.Add(new int[] {4,14, 21});
                    CS.Add(new int[] {4,14, 20});

                    CS.Add(new int[] {3,15, 22});
                    CS.Add(new int[] {3,15, 21});
                    CS.Add(new int[] {3,15, 20});

                    CS.Add(new int[] {4,15, 22});
                    CS.Add(new int[] {4,15, 21});
                    CS.Add(new int[] {4,15, 20});

                    CS.Add(new int[] {7,15, 13,19});
                    CS.Add(new int[] {6,15, 13,19});
                    CS.Add(new int[] {5,15, 13,19});

                    CS.Add(new int[] {7,14, 13,19});
                    CS.Add(new int[] {6,14, 13,19});
                    CS.Add(new int[] {5,14, 13,19});

                    CS.Add(new int[] {7,15, 12,19});
                    CS.Add(new int[] {6,15, 12,19});
                    CS.Add(new int[] {5,15, 12,19});

                    CS.Add(new int[] {7,15, 11,19});
                    CS.Add(new int[] {6,15, 11,19});
                    CS.Add(new int[] {5,15, 11,19});

                    CS.Add(new int[] {7,15, 10,19});
                    CS.Add(new int[] {6,15, 10,19});
                    CS.Add(new int[] {5,15, 10,19});



                    CS.Add(new int[] {7,14, 12,19});
                    CS.Add(new int[] {6,14, 12,19});
                    CS.Add(new int[] {5,14, 12,19});

                    CS.Add(new int[] {7,14, 11,19});
                    CS.Add(new int[] {6,14, 11,19});
                    CS.Add(new int[] {5,14, 11,19});

                    CS.Add(new int[] {7,14, 10,19});
                    CS.Add(new int[] {6,14, 10,19});
                    CS.Add(new int[] {5,14, 10,19});

                    CS.Add(new int[] {7,15, 13,18});
                    CS.Add(new int[] {6,15, 13,18});
                    CS.Add(new int[] {5,15, 13,18});

                    CS.Add(new int[] {7,14, 13,18});
                    CS.Add(new int[] {6,14, 13,18});
                    CS.Add(new int[] {5,14, 13,18});

                    CS.Add(new int[] {7,15, 12,18});
                    CS.Add(new int[] {6,15, 12,18});
                    CS.Add(new int[] {5,15, 12,18});

                    CS.Add(new int[] {7,15, 11,18});
                    CS.Add(new int[] {6,15, 11,18});
                    CS.Add(new int[] {5,15, 11,18});

                    CS.Add(new int[] {7,15, 10,18});
                    CS.Add(new int[] {6,15, 10,18});
                    CS.Add(new int[] {5,15, 10,18});

                    CS.Add(new int[] {7,14, 12,19});
                    CS.Add(new int[] {6,14, 12,19});
                    CS.Add(new int[] {5,14, 12,19});

                    CS.Add(new int[] {7,14, 11,18});
                    CS.Add(new int[] {6,14, 11,18});
                    CS.Add(new int[] {5,14, 11,18});

                    CS.Add(new int[] {7,14, 10,18});
                    CS.Add(new int[] {6,14, 10,18});
                    CS.Add(new int[] {5,14, 10,18});

                    CS.Add(new int[] {3, 8,10, 16});
                    CS.Add(new int[] {3, 8,11, 16});
                    CS.Add(new int[] {3, 8,12, 16});
                    CS.Add(new int[] {3, 8,13, 16});

                    CS.Add(new int[] {3, 9,10, 16});
                    CS.Add(new int[] {3, 9,11, 16});
                    CS.Add(new int[] {3, 9,12, 16});
                    CS.Add(new int[] {3, 9,13, 16});

                    CS.Add(new int[] {4, 8,10, 16});
                    CS.Add(new int[] {4, 8,11, 16});
                    CS.Add(new int[] {4, 8,12, 16});
                    CS.Add(new int[] {4, 8,13, 16});

                    CS.Add(new int[] {4, 9,10, 16});
                    CS.Add(new int[] {4, 9,11, 16});
                    CS.Add(new int[] {4, 9,12, 16});
                    CS.Add(new int[] {4, 9,13, 16});


                    CS.Add(new int[] {3, 8,10, 17});
                    CS.Add(new int[] {3, 8,11, 17});
                    CS.Add(new int[] {3, 8,12, 17});
                    CS.Add(new int[] {3, 8,13, 17});

                    CS.Add(new int[] {3, 9,10, 17});
                    CS.Add(new int[] {3, 9,11, 17});
                    CS.Add(new int[] {3, 9,12, 17});
                    CS.Add(new int[] {3, 9,13, 17});

                    CS.Add(new int[] {4, 8,10, 17});
                    CS.Add(new int[] {4, 8,11, 17});
                    CS.Add(new int[] {4, 8,12, 17});
                    CS.Add(new int[] {4, 8,13, 17});

                    CS.Add(new int[] {4, 9,10, 17});
                    CS.Add(new int[] {4, 9,11, 17});
                    CS.Add(new int[] {4, 9,12, 17});
                    CS.Add(new int[] {4, 9,13, 17});

                    CS.Add(new int[] {3, 14, 20});
                    CS.Add(new int[] {3, 14, 21});
                    CS.Add(new int[] {3, 14, 22});

                    CS.Add(new int[] {3, 15, 20});
                    CS.Add(new int[] {3, 15, 21});
                    CS.Add(new int[] {3, 15, 22});

                    CS.Add(new int[] {4, 14, 20});
                    CS.Add(new int[] {4, 14, 21});
                    CS.Add(new int[] {4, 14, 22});

                    CS.Add(new int[] {4, 15, 20});
                    CS.Add(new int[] {4, 15, 21});
                    CS.Add(new int[] {4, 15, 22});

                    CS.Add(new int[] {3, 10, 18});
                    CS.Add(new int[] {3, 11, 18});
                    CS.Add(new int[] {3, 12, 18});
                    CS.Add(new int[] {3, 13, 18});

                    CS.Add(new int[] {4, 10, 18});
                    CS.Add(new int[] {4, 11, 18});
                    CS.Add(new int[] {4, 12, 18});
                    CS.Add(new int[] {4, 13, 18});

                    CS.Add(new int[] {3, 10, 19});
                    CS.Add(new int[] {3, 11, 19});
                    CS.Add(new int[] {3, 12, 19});
                    CS.Add(new int[] {3, 13, 19});

                    CS.Add(new int[] {4, 10, 19});
                    CS.Add(new int[] {4, 11, 19});
                    CS.Add(new int[] {4, 12, 19});
                    CS.Add(new int[] {4, 13, 19});

                    PS.Add(new int[] {1, 2, 3, 4, 5, 6, 7});
                    PS.Add(new int[] {1, 2, 3, 4, 14, 15, 20, 21, 22});
                    PS.Add(new int[] {1, 2, 5, 6, 7, 10, 11, 12, 13, 14, 15});
                    PS.Add(new int[] {1, 2, 10, 11, 12, 13, 20, 21, 22});
                    PS.Add(new int[] {1, 2, 8, 9, 18, 19, 20, 21, 22});
                    PS.Add(new int[] {1, 2, 5, 6, 7, 8, 9, 14, 15, 18, 19});
                    PS.Add(new int[] {3, 4, 5, 6, 7, 8, 9, 16, 17});
                    PS.Add(new int[] {3, 4, 8, 9, 14, 15, 16, 17, 20, 21, 22});
                    PS.Add(new int[] {8, 9, 10, 11, 12, 13, 16, 17, 20, 21, 22});
                    PS.Add(new int[] {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17});
                    PS.Add(new int[] {3, 4, 5, 6, 7, 10, 11, 12, 13, 16, 17, 18, 19});
                    PS.Add(new int[] {5, 6, 7, 14, 15, 16, 17, 18, 19});
                    PS.Add(new int[] {16, 17, 18, 19, 20, 21, 22});
                    break;
                case "2o3":
                    number = 3;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7};
                    /*correlations = new double[] { 1.00, -0.28, 0.1,
                                                -0.28,  1.00, 0.2,
                                                0.1,   0.2,  1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228,
                                                  0.1377770, 1.0000000, 0.1108590,
                                                  0.0326228, 0.1108590, 1.0000000};
                    correlations_matrix = Matrix<double>.Build.Dense(3,3,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2});
                    CS.Add(new int[] {1,3});
                    CS.Add(new int[] {2,3});
                    PS.Add(new int[] {1,2});
                    PS.Add(new int[] {1,3});
                    PS.Add(new int[] {2,3});
                    break;
                case "3o5":
                    number = 5;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.0719920, -0.0486818,
                                                  0.1377770, 1.0000000, 0.1108590, 0.0981194,  0.1800180,
                                                  0.0326228, 0.1108590, 1.0000000, 0.0873844,  0.1303380,
                                                  0.0719920, 0.0981194, 0.0873844, 1.0000000,  0.0248845, 
                                                 -0.0486818, 0.1800180, 0.1303380, 0.0248845,  1.0000000};  

                    correlations_matrix = Matrix<double>.Build.Dense(5,5,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2,3});
                    CS.Add(new int[] {1,2,4});
                    CS.Add(new int[] {1,2,5});
                    CS.Add(new int[] {1,3,4});
                    CS.Add(new int[] {1,3,5});
                    CS.Add(new int[] {1,4,5});
                    CS.Add(new int[] {2,3,4});
                    CS.Add(new int[] {2,3,5});
                    CS.Add(new int[] {2,4,5});
                    CS.Add(new int[] {3,4,5});

                    PS.Add(new int[] {1,2,3});
                    PS.Add(new int[] {1,2,4});
                    PS.Add(new int[] {1,2,5});
                    PS.Add(new int[] {1,3,4});
                    PS.Add(new int[] {1,3,5});
                    PS.Add(new int[] {1,4,5});
                    PS.Add(new int[] {2,3,4});
                    PS.Add(new int[] {2,3,5});
                    PS.Add(new int[] {2,4,5});
                    PS.Add(new int[] {3,4,5});
                    break;
                case "4o7":
                    number = 7;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77, 0.67};
                    /*correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00};*/
                    /*correlations = new double[] { 1.00, 0, 0, 0, 0, 0, 0,
                                                    0, 1.00, 0, 0, 0, 0, 0,
                                                    0, 0, 1.00, 0, 0, 0, 0,
                                                    0, 0, 0, 1.00, 0, 0, 0, 
                                                    0, 0, 0, 0, 1.00, 0, 0,
                                                    0, 0, 0, 0, 0, 1.00, 0,
                                                    0, 0, 0, 0, 0, 0, 1.00};*/
                    correlations = new double[] { 1.0000000, 0.1377770, 0.0326228, 0.071992000, -0.04868180, 0.09131980, -0.066384200, 
                                                  0.1377770, 1.0000000, 0.1108590, 0.098119400,  0.18001800, 0.01634070, -0.059739400, 
                                                  0.0326228, 0.1108590, 1.0000000, 0.087384400,  0.13033800, 0.11618600, -0.056192500, 
                                                  0.0719920, 0.0981194, 0.0873844, 1.000000000,  0.02488450,-0.03212350, -0.000253982, 
                                                 -0.0486818, 0.1800180, 0.1303380, 0.024884500,  1.00000000,-0.00696747,  0.127017000, 
                                                  0.0913198, 0.0163407, 0.1161860,-0.032123500, -0.00696747, 1.00000000,  0.014421100, 
                                                 -0.0663842,-0.0597394,-0.0561925,-0.000253982,  0.12701700, 0.01442110,  1.000000000};
                    correlations_matrix = Matrix<double>.Build.Dense(7,7,correlations);
                    GenerateSigma();
                    List<int[]> wolfram4 = new List<int[]>();
                    wolfram4.AddRange(new List<int[]> {new int[] {1, 1, 1, 1, 0, 0, 0}, 
                                    new int[] {1, 1, 1, 0, 1, 0, 0}, 
                                    new int[] {1, 1, 1, 0, 0, 1, 0}, 
                                    new int[] {1, 1, 1, 0, 0, 0, 1}, 
                                    new int[] {1, 1, 0, 1, 1, 0, 0}, 
                                    new int[] {1, 1, 0, 1, 0, 1, 0}, 
                                    new int[] {1, 1, 0, 1, 0, 0, 1}, 
                                    new int[] {1, 1, 0, 0, 1, 1, 0}, 
                                    new int[] {1, 1, 0, 0, 1, 0, 1}, 
                                    new int[] {1, 1, 0, 0, 0, 1, 1}, 
                                    new int[] {1, 0, 1, 1, 1, 0, 0}, 
                                    new int[] {1, 0, 1, 1, 0, 1, 0}, 
                                    new int[] {1, 0, 1, 1, 0, 0, 1}, 
                                    new int[] {1, 0, 1, 0, 1, 1, 0}, 
                                    new int[] {1, 0, 1, 0, 1, 0, 1}, 
                                    new int[] {1, 0, 1, 0, 0, 1, 1}, 
                                    new int[] {1, 0, 0, 1, 1, 1, 0}, 
                                    new int[] {1, 0, 0, 1, 1, 0, 1}, 
                                    new int[] {1, 0, 0, 1, 0, 1, 1}, 
                                    new int[] {1, 0, 0, 0, 1, 1, 1}, 
                                    new int[] {0, 1, 1, 1, 1, 0, 0}, 
                                    new int[] {0, 1, 1, 1, 0, 1, 0}, 
                                    new int[] {0, 1, 1, 1, 0, 0, 1}, 
                                    new int[] {0, 1, 1, 0, 1, 1, 0}, 
                                    new int[] {0, 1, 1, 0, 1, 0, 1}, 
                                    new int[] {0, 1, 1, 0, 0, 1, 1}, 
                                    new int[] {0, 1, 0, 1, 1, 1, 0}, 
                                    new int[] {0, 1, 0, 1, 1, 0, 1}, 
                                    new int[] {0, 1, 0, 1, 0, 1, 1}, 
                                    new int[] {0, 1, 0, 0, 1, 1, 1}, 
                                    new int[] {0, 0, 1, 1, 1, 1, 0}, 
                                    new int[] {0, 0, 1, 1, 1, 0, 1}, 
                                    new int[] {0, 0, 1, 1, 0, 1, 1}, 
                                    new int[] {0, 0, 1, 0, 1, 1, 1}, 
                                    new int[] {0, 0, 0, 1, 1, 1, 1}});
                    CS = new List<int[]>(WolframPerms(wolfram4));
                    PS = new List<int[]>(CS);
                    break;
                case "6o10":
                    number = 10;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7, 0.72, 0.73, 0.77, 0.67, 0.96, 0.72, 0.88};
                    correlations = new double[] {   1.00, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 1.00, 0, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 0, 1.00, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 1.00, 0, 0, 0, 0, 0, 0, 
                                                    0, 0, 0, 0, 1.00, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 1.00, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0, 1.00, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0, 0, 1.00, 0, 0,
                                                    0, 0, 0, 0, 0, 0, 0, 0, 1.00, 0,
                                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1.00};
                    correlations_matrix = Matrix<double>.Build.Dense(10,10,correlations);
                    GenerateSigma();
                    List<int[]> wolfram6 = new List<int[]>();
                    wolfram6.AddRange(new List<int[]> { new int[] {1, 1, 1, 1, 1, 1, 0, 0, 0, 0}, 
                                                        new int[] {1, 1, 1, 1, 1, 0, 1, 0, 0, 0}, 
                                                        new int[] {1, 1, 1, 1, 1, 0, 0, 1, 0, 0}, 
                                                        new int[] {1, 1, 1, 1, 1, 0, 0, 0, 1, 0}, 
                                                        new int[] {1, 1, 1, 1, 1, 0, 0, 0, 0, 1}, 
                                                        new int[] {1, 1, 1, 1, 0, 1, 1, 0, 0, 0}, 
                                                        new int[] {1, 1, 1, 1, 0, 1, 0, 1, 0, 0}, 
                                                        new int[] {1, 1, 1, 1, 0, 1, 0, 0, 1, 0}, 
                                                        new int[] {1, 1, 1, 1, 0, 1, 0, 0, 0, 1}, 
                                                        new int[] {1, 1, 1, 1, 0, 0, 1, 1, 0, 0}, 
                                                        new int[] {1, 1, 1, 1, 0, 0, 1, 0, 1, 0}, 
                                                        new int[] {1, 1, 1, 1, 0, 0, 1, 0, 0, 1}, 
                                                        new int[] {1, 1, 1, 1, 0, 0, 0, 1, 1, 0}, 
                                                        new int[] {1, 1, 1, 1, 0, 0, 0, 1, 0, 1}, 
                                                        new int[] {1, 1, 1, 1, 0, 0, 0, 0, 1, 1}, 
                                                        new int[] {1, 1, 1, 0, 1, 1, 1, 0, 0, 0}, 
                                                        new int[] {1, 1, 1, 0, 1, 1, 0, 1, 0, 0}, 
                                                        new int[] {1, 1, 1, 0, 1, 1, 0, 0, 1, 0}, 
                                                        new int[] {1, 1, 1, 0, 1, 1, 0, 0, 0, 1}, 
                                                        new int[] {1, 1, 1, 0, 1, 0, 1, 1, 0, 0}, 
                                                        new int[] {1, 1, 1, 0, 1, 0, 1, 0, 1, 0}, 
                                                        new int[] {1, 1, 1, 0, 1, 0, 1, 0, 0, 1}, 
                                                        new int[] {1, 1, 1, 0, 1, 0, 0, 1, 1, 0}, 
                                                        new int[] {1, 1, 1, 0, 1, 0, 0, 1, 0, 1}, 
                                                        new int[] {1, 1, 1, 0, 1, 0, 0, 0, 1, 1}, 
                                                        new int[] {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, 
                                                        new int[] {1, 1, 1, 0, 0, 1, 1, 0, 1, 0}, 
                                                        new int[] {1, 1, 1, 0, 0, 1, 1, 0, 0, 1}, 
                                                        new int[] {1, 1, 1, 0, 0, 1, 0, 1, 1, 0}, 
                                                        new int[] {1, 1, 1, 0, 0, 1, 0, 1, 0, 1}, 
                                                        new int[] {1, 1, 1, 0, 0, 1, 0, 0, 1, 1}, 
                                                        new int[] {1, 1, 1, 0, 0, 0, 1, 1, 1, 0}, 
                                                        new int[] {1, 1, 1, 0, 0, 0, 1, 1, 0, 1}, 
                                                        new int[] {1, 1, 1, 0, 0, 0, 1, 0, 1, 1}, 
                                                        new int[] {1, 1, 1, 0, 0, 0, 0, 1, 1, 1}, 
                                                        new int[] {1, 1, 0, 1, 1, 1, 1, 0, 0, 0}, 
                                                        new int[] {1, 1, 0, 1, 1, 1, 0, 1, 0, 0}, 
                                                        new int[] {1, 1, 0, 1, 1, 1, 0, 0, 1, 0}, 
                                                        new int[] {1, 1, 0, 1, 1, 1, 0, 0, 0, 1}, 
                                                        new int[] {1, 1, 0, 1, 1, 0, 1, 1, 0, 0}, 
                                                        new int[] {1, 1, 0, 1, 1, 0, 1, 0, 1, 0}, 
                                                        new int[] {1, 1, 0, 1, 1, 0, 1, 0, 0, 1}, 
                                                        new int[] {1, 1, 0, 1, 1, 0, 0, 1, 1, 0}, 
                                                        new int[] {1, 1, 0, 1, 1, 0, 0, 1, 0, 1}, 
                                                        new int[] {1, 1, 0, 1, 1, 0, 0, 0, 1, 1}, 
                                                        new int[] {1, 1, 0, 1, 0, 1, 1, 1, 0, 0}, 
                                                        new int[] {1, 1, 0, 1, 0, 1, 1, 0, 1, 0}, 
                                                        new int[] {1, 1, 0, 1, 0, 1, 1, 0, 0, 1}, 
                                                        new int[] {1, 1, 0, 1, 0, 1, 0, 1, 1, 0}, 
                                                        new int[] {1, 1, 0, 1, 0, 1, 0, 1, 0, 1}, 
                                                        new int[] {1, 1, 0, 1, 0, 1, 0, 0, 1, 1}, 
                                                        new int[] {1, 1, 0, 1, 0, 0, 1, 1, 1, 0}, 
                                                        new int[] {1, 1, 0, 1, 0, 0, 1, 1, 0, 1}, 
                                                        new int[] {1, 1, 0, 1, 0, 0, 1, 0, 1, 1}, 
                                                        new int[] {1, 1, 0, 1, 0, 0, 0, 1, 1, 1}, 
                                                        new int[] {1, 1, 0, 0, 1, 1, 1, 1, 0, 0}, 
                                                        new int[] {1, 1, 0, 0, 1, 1, 1, 0, 1, 0}, 
                                                        new int[] {1, 1, 0, 0, 1, 1, 1, 0, 0, 1}, 
                                                        new int[] {1, 1, 0, 0, 1, 1, 0, 1, 1, 0}, 
                                                        new int[] {1, 1, 0, 0, 1, 1, 0, 1, 0, 1}, 
                                                        new int[] {1, 1, 0, 0, 1, 1, 0, 0, 1, 1}, 
                                                        new int[] {1, 1, 0, 0, 1, 0, 1, 1, 1, 0}, 
                                                        new int[] {1, 1, 0, 0, 1, 0, 1, 1, 0, 1}, 
                                                        new int[] {1, 1, 0, 0, 1, 0, 1, 0, 1, 1}, 
                                                        new int[] {1, 1, 0, 0, 1, 0, 0, 1, 1, 1}, 
                                                        new int[] {1, 1, 0, 0, 0, 1, 1, 1, 1, 0}, 
                                                        new int[] {1, 1, 0, 0, 0, 1, 1, 1, 0, 1}, 
                                                        new int[] {1, 1, 0, 0, 0, 1, 1, 0, 1, 1}, 
                                                        new int[] {1, 1, 0, 0, 0, 1, 0, 1, 1, 1}, 
                                                        new int[] {1, 1, 0, 0, 0, 0, 1, 1, 1, 1}, 
                                                        new int[] {1, 0, 1, 1, 1, 1, 1, 0, 0, 0}, 
                                                        new int[] {1, 0, 1, 1, 1, 1, 0, 1, 0, 0}, 
                                                        new int[] {1, 0, 1, 1, 1, 1, 0, 0, 1, 0}, 
                                                        new int[] {1, 0, 1, 1, 1, 1, 0, 0, 0, 1}, 
                                                        new int[] {1, 0, 1, 1, 1, 0, 1, 1, 0, 0}, 
                                                        new int[] {1, 0, 1, 1, 1, 0, 1, 0, 1, 0}, 
                                                        new int[] {1, 0, 1, 1, 1, 0, 1, 0, 0, 1}, 
                                                        new int[] {1, 0, 1, 1, 1, 0, 0, 1, 1, 0}, 
                                                        new int[] {1, 0, 1, 1, 1, 0, 0, 1, 0, 1}, 
                                                        new int[] {1, 0, 1, 1, 1, 0, 0, 0, 1, 1}, 
                                                        new int[] {1, 0, 1, 1, 0, 1, 1, 1, 0, 0}, 
                                                        new int[] {1, 0, 1, 1, 0, 1, 1, 0, 1, 0}, 
                                                        new int[] {1, 0, 1, 1, 0, 1, 1, 0, 0, 1}, 
                                                        new int[] {1, 0, 1, 1, 0, 1, 0, 1, 1, 0}, 
                                                        new int[] {1, 0, 1, 1, 0, 1, 0, 1, 0, 1}, 
                                                        new int[] {1, 0, 1, 1, 0, 1, 0, 0, 1, 1}, 
                                                        new int[] {1, 0, 1, 1, 0, 0, 1, 1, 1, 0}, 
                                                        new int[] {1, 0, 1, 1, 0, 0, 1, 1, 0, 1}, 
                                                        new int[] {1, 0, 1, 1, 0, 0, 1, 0, 1, 1}, 
                                                        new int[] {1, 0, 1, 1, 0, 0, 0, 1, 1, 1}, 
                                                        new int[] {1, 0, 1, 0, 1, 1, 1, 1, 0, 0}, 
                                                        new int[] {1, 0, 1, 0, 1, 1, 1, 0, 1, 0}, 
                                                        new int[] {1, 0, 1, 0, 1, 1, 1, 0, 0, 1}, 
                                                        new int[] {1, 0, 1, 0, 1, 1, 0, 1, 1, 0}, 
                                                        new int[] {1, 0, 1, 0, 1, 1, 0, 1, 0, 1}, 
                                                        new int[] {1, 0, 1, 0, 1, 1, 0, 0, 1, 1}, 
                                                        new int[] {1, 0, 1, 0, 1, 0, 1, 1, 1, 0}, 
                                                        new int[] {1, 0, 1, 0, 1, 0, 1, 1, 0, 1}, 
                                                        new int[] {1, 0, 1, 0, 1, 0, 1, 0, 1, 1}, 
                                                        new int[] {1, 0, 1, 0, 1, 0, 0, 1, 1, 1}, 
                                                        new int[] {1, 0, 1, 0, 0, 1, 1, 1, 1, 0}, 
                                                        new int[] {1, 0, 1, 0, 0, 1, 1, 1, 0, 1}, 
                                                        new int[] {1, 0, 1, 0, 0, 1, 1, 0, 1, 1}, 
                                                        new int[] {1, 0, 1, 0, 0, 1, 0, 1, 1, 1}, 
                                                        new int[] {1, 0, 1, 0, 0, 0, 1, 1, 1, 1}, 
                                                        new int[] {1, 0, 0, 1, 1, 1, 1, 1, 0, 0}, 
                                                        new int[] {1, 0, 0, 1, 1, 1, 1, 0, 1, 0}, 
                                                        new int[] {1, 0, 0, 1, 1, 1, 1, 0, 0, 1}, 
                                                        new int[] {1, 0, 0, 1, 1, 1, 0, 1, 1, 0}, 
                                                        new int[] {1, 0, 0, 1, 1, 1, 0, 1, 0, 1}, 
                                                        new int[] {1, 0, 0, 1, 1, 1, 0, 0, 1, 1}, 
                                                        new int[] {1, 0, 0, 1, 1, 0, 1, 1, 1, 0}, 
                                                        new int[] {1, 0, 0, 1, 1, 0, 1, 1, 0, 1}, 
                                                        new int[] {1, 0, 0, 1, 1, 0, 1, 0, 1, 1}, 
                                                        new int[] {1, 0, 0, 1, 1, 0, 0, 1, 1, 1}, 
                                                        new int[] {1, 0, 0, 1, 0, 1, 1, 1, 1, 0}, 
                                                        new int[] {1, 0, 0, 1, 0, 1, 1, 1, 0, 1}, 
                                                        new int[] {1, 0, 0, 1, 0, 1, 1, 0, 1, 1}, 
                                                        new int[] {1, 0, 0, 1, 0, 1, 0, 1, 1, 1}, 
                                                        new int[] {1, 0, 0, 1, 0, 0, 1, 1, 1, 1}, 
                                                        new int[] {1, 0, 0, 0, 1, 1, 1, 1, 1, 0}, 
                                                        new int[] {1, 0, 0, 0, 1, 1, 1, 1, 0, 1}, 
                                                        new int[] {1, 0, 0, 0, 1, 1, 1, 0, 1, 1}, 
                                                        new int[] {1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 
                                                        new int[] {1, 0, 0, 0, 1, 0, 1, 1, 1, 1}, 
                                                        new int[] {1, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
                                                        new int[] {0, 1, 1, 1, 1, 1, 1, 0, 0, 0}, 
                                                        new int[] {0, 1, 1, 1, 1, 1, 0, 1, 0, 0}, 
                                                        new int[] {0, 1, 1, 1, 1, 1, 0, 0, 1, 0}, 
                                                        new int[] {0, 1, 1, 1, 1, 1, 0, 0, 0, 1}, 
                                                        new int[] {0, 1, 1, 1, 1, 0, 1, 1, 0, 0}, 
                                                        new int[] {0, 1, 1, 1, 1, 0, 1, 0, 1, 0}, 
                                                        new int[] {0, 1, 1, 1, 1, 0, 1, 0, 0, 1}, 
                                                        new int[] {0, 1, 1, 1, 1, 0, 0, 1, 1, 0}, 
                                                        new int[] {0, 1, 1, 1, 1, 0, 0, 1, 0, 1}, 
                                                        new int[] {0, 1, 1, 1, 1, 0, 0, 0, 1, 1}, 
                                                        new int[] {0, 1, 1, 1, 0, 1, 1, 1, 0, 0}, 
                                                        new int[] {0, 1, 1, 1, 0, 1, 1, 0, 1, 0}, 
                                                        new int[] {0, 1, 1, 1, 0, 1, 1, 0, 0, 1}, 
                                                        new int[] {0, 1, 1, 1, 0, 1, 0, 1, 1, 0}, 
                                                        new int[] {0, 1, 1, 1, 0, 1, 0, 1, 0, 1}, 
                                                        new int[] {0, 1, 1, 1, 0, 1, 0, 0, 1, 1}, 
                                                        new int[] {0, 1, 1, 1, 0, 0, 1, 1, 1, 0}, 
                                                        new int[] {0, 1, 1, 1, 0, 0, 1, 1, 0, 1}, 
                                                        new int[] {0, 1, 1, 1, 0, 0, 1, 0, 1, 1}, 
                                                        new int[] {0, 1, 1, 1, 0, 0, 0, 1, 1, 1}, 
                                                        new int[] {0, 1, 1, 0, 1, 1, 1, 1, 0, 0}, 
                                                        new int[] {0, 1, 1, 0, 1, 1, 1, 0, 1, 0}, 
                                                        new int[] {0, 1, 1, 0, 1, 1, 1, 0, 0, 1}, 
                                                        new int[] {0, 1, 1, 0, 1, 1, 0, 1, 1, 0}, 
                                                        new int[] {0, 1, 1, 0, 1, 1, 0, 1, 0, 1}, 
                                                        new int[] {0, 1, 1, 0, 1, 1, 0, 0, 1, 1}, 
                                                        new int[] {0, 1, 1, 0, 1, 0, 1, 1, 1, 0}, 
                                                        new int[] {0, 1, 1, 0, 1, 0, 1, 1, 0, 1}, 
                                                        new int[] {0, 1, 1, 0, 1, 0, 1, 0, 1, 1}, 
                                                        new int[] {0, 1, 1, 0, 1, 0, 0, 1, 1, 1}, 
                                                        new int[] {0, 1, 1, 0, 0, 1, 1, 1, 1, 0}, 
                                                        new int[] {0, 1, 1, 0, 0, 1, 1, 1, 0, 1}, 
                                                        new int[] {0, 1, 1, 0, 0, 1, 1, 0, 1, 1}, 
                                                        new int[] {0, 1, 1, 0, 0, 1, 0, 1, 1, 1}, 
                                                        new int[] {0, 1, 1, 0, 0, 0, 1, 1, 1, 1}, 
                                                        new int[] {0, 1, 0, 1, 1, 1, 1, 1, 0, 0}, 
                                                        new int[] {0, 1, 0, 1, 1, 1, 1, 0, 1, 0}, 
                                                        new int[] {0, 1, 0, 1, 1, 1, 1, 0, 0, 1}, 
                                                        new int[] {0, 1, 0, 1, 1, 1, 0, 1, 1, 0}, 
                                                        new int[] {0, 1, 0, 1, 1, 1, 0, 1, 0, 1}, 
                                                        new int[] {0, 1, 0, 1, 1, 1, 0, 0, 1, 1}, 
                                                        new int[] {0, 1, 0, 1, 1, 0, 1, 1, 1, 0}, 
                                                        new int[] {0, 1, 0, 1, 1, 0, 1, 1, 0, 1}, 
                                                        new int[] {0, 1, 0, 1, 1, 0, 1, 0, 1, 1}, 
                                                        new int[] {0, 1, 0, 1, 1, 0, 0, 1, 1, 1}, 
                                                        new int[] {0, 1, 0, 1, 0, 1, 1, 1, 1, 0}, 
                                                        new int[] {0, 1, 0, 1, 0, 1, 1, 1, 0, 1}, 
                                                        new int[] {0, 1, 0, 1, 0, 1, 1, 0, 1, 1}, 
                                                        new int[] {0, 1, 0, 1, 0, 1, 0, 1, 1, 1}, 
                                                        new int[] {0, 1, 0, 1, 0, 0, 1, 1, 1, 1}, 
                                                        new int[] {0, 1, 0, 0, 1, 1, 1, 1, 1, 0}, 
                                                        new int[] {0, 1, 0, 0, 1, 1, 1, 1, 0, 1}, 
                                                        new int[] {0, 1, 0, 0, 1, 1, 1, 0, 1, 1}, 
                                                        new int[] {0, 1, 0, 0, 1, 1, 0, 1, 1, 1}, 
                                                        new int[] {0, 1, 0, 0, 1, 0, 1, 1, 1, 1}, 
                                                        new int[] {0, 1, 0, 0, 0, 1, 1, 1, 1, 1}, 
                                                        new int[] {0, 0, 1, 1, 1, 1, 1, 1, 0, 0}, 
                                                        new int[] {0, 0, 1, 1, 1, 1, 1, 0, 1, 0}, 
                                                        new int[] {0, 0, 1, 1, 1, 1, 1, 0, 0, 1}, 
                                                        new int[] {0, 0, 1, 1, 1, 1, 0, 1, 1, 0}, 
                                                        new int[] {0, 0, 1, 1, 1, 1, 0, 1, 0, 1}, 
                                                        new int[] {0, 0, 1, 1, 1, 1, 0, 0, 1, 1}, 
                                                        new int[] {0, 0, 1, 1, 1, 0, 1, 1, 1, 0}, 
                                                        new int[] {0, 0, 1, 1, 1, 0, 1, 1, 0, 1}, 
                                                        new int[] {0, 0, 1, 1, 1, 0, 1, 0, 1, 1}, 
                                                        new int[] {0, 0, 1, 1, 1, 0, 0, 1, 1, 1}, 
                                                        new int[] {0, 0, 1, 1, 0, 1, 1, 1, 1, 0}, 
                                                        new int[] {0, 0, 1, 1, 0, 1, 1, 1, 0, 1}, 
                                                        new int[] {0, 0, 1, 1, 0, 1, 1, 0, 1, 1}, 
                                                        new int[] {0, 0, 1, 1, 0, 1, 0, 1, 1, 1}, 
                                                        new int[] {0, 0, 1, 1, 0, 0, 1, 1, 1, 1}, 
                                                        new int[] {0, 0, 1, 0, 1, 1, 1, 1, 1, 0}, 
                                                        new int[] {0, 0, 1, 0, 1, 1, 1, 1, 0, 1}, 
                                                        new int[] {0, 0, 1, 0, 1, 1, 1, 0, 1, 1}, 
                                                        new int[] {0, 0, 1, 0, 1, 1, 0, 1, 1, 1}, 
                                                        new int[] {0, 0, 1, 0, 1, 0, 1, 1, 1, 1}, 
                                                        new int[] {0, 0, 1, 0, 0, 1, 1, 1, 1, 1}, 
                                                        new int[] {0, 0, 0, 1, 1, 1, 1, 1, 1, 0}, 
                                                        new int[] {0, 0, 0, 1, 1, 1, 1, 1, 0, 1}, 
                                                        new int[] {0, 0, 0, 1, 1, 1, 1, 0, 1, 1}, 
                                                        new int[] {0, 0, 0, 1, 1, 1, 0, 1, 1, 1}, 
                                                        new int[] {0, 0, 0, 1, 1, 0, 1, 1, 1, 1}, 
                                                        new int[] {0, 0, 0, 1, 0, 1, 1, 1, 1, 1}, 
                                                        new int[] {0, 0, 0, 0, 1, 1, 1, 1, 1, 1}});
                    CS = new List<int[]>(WolframPerms(wolfram6));
                    PS = new List<int[]>(CS);
                    break;
                default:
                    break;
            }
        }
       
        static void Generate_Bm(int v = 0)
        {
            List<double> temp = new List<double>();
            if (v == 0)
            {
                for(int i = 2; i <= number; i++)
                {
                    temp.Clear();
                    temp.Add(0);
                    temp.AddRange(b(i));
                    Bm.Add(temp.ToArray());
                }
            }
            else
            {
                    temp.Clear();
                    temp.Add(0);
                    temp.AddRange(b(v));
                    Bm.Add(temp.ToArray());
            }

        }

        //Needs further optimization with hashing at node values to eliminate redundant computations
        static double Simulation(int runs)
        {
            List<int> Path = new List<int>();           //Path used while simulating tree traversal
            List<int[]> Paths = RecursiveEnumerate();   //List of Paths (used to determine cutsets)
            Random Rand = new Random();                 //Random used to determine left or right traversal
            Hashtable ht = new Hashtable();             //Hashtable used to lookup node prob already calculated before
            List<PTrace> sim = new List<PTrace>();      //Simulation Graph
            List<PTrace> num = new List<PTrace>();      //Numerical Graph
            List<PTrace> err = new List<PTrace>();      //Difference between Simulation and Numerical Graph
            Generate_Bm();
            //var pkey = Path.Concat(new List<int> {1}).ToArray(); //path key used in hashtable
            double choice;
            double condi;
            int success = 0;

            for (int i = 0; i < runs; i++)
            {
                Path.Clear();
                for (int j = 0; j < number; j++)
                {
                    choice = Rand.NextDouble();
                    //If this is the first component in the path
                    if (j == 0)
                    {
                        if (choice < component_reliabilities[0])
                            Path.Add(1);
                        else
                            Path.Add(0);

                        if (!sim.Contains(sim.Find(pT => pT.path == String.Join("",Path.ToArray()))))
                            sim.Add(new PTrace(String.Join("",Path.ToArray())));
                        else
                            sim.Find(pT => pT.path == String.Join("",Path.ToArray())).count++;
                    }
                    else
                    {
                        //Check if we have the B vector for the upcoming calculation
                        /*if (Bm.Count < Path.Count + 2)  //Bm = {0,0,...}
                            Generate_Bm(Path.Count + 1);//Bm = {0,0,x,...}*/

                        //Check if prob is in hash
                        if (!ht.ContainsKey(String.Join("",Path.Concat(new List<int> {1}).ToArray())))
                        {
                            condi = ConditionalProbability(Path.Concat(new List<int> {1}).ToArray());
                            ht.Add(String.Join("",Path.Concat(new List<int> {1}).ToArray()), condi);
                        }

                        if (choice < (double)ht[String.Join("",Path.Concat(new List<int> {1}).ToArray())])
                            Path.Add(1);
                        else
                            Path.Add(0);

                        if (!sim.Contains(sim.Find(pT => pT.path == String.Join("",Path.ToArray()))))
                            sim.Add(new PTrace(String.Join("",Path.ToArray())));
                        else
                            sim.Find(pT => pT.path == String.Join("",Path.ToArray())).count++;
                    }
                    if (ContainsCutset(Path.ToArray()))
                        break;
                    if (ContainsPathset(Path.ToArray()))
                    {
                        success++;
                        break;
                    }
                }
                
            }
            /* //Create Graphviz representations

            foreach (PTrace p in sim)
            {
                p.prob = (double) p.count / (double) runs;
                num.Add(new PTrace(p.path, TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c)))));
                err.Add(new PTrace(p.path, Math.Abs(p.prob - TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c))))));
                
                //if (Math.Abs(p.prob - TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c)))) > 0.001)
                   // totals.Add(new PTrace(p.path, Math.Abs(p.prob - TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c))))));
            }
            CreateDOTGraphs(sim,num,err); */

            return (double) success / (double) runs; 
        }

        static void CreateDOTGraphs(List<PTrace> s, List<PTrace> n, List<PTrace> e)
        {
            string sim_path = @"Z:\simulation.dot";
            string num_path = @"Z:\numerical.dot";
            string err_path = @"Z:\error.dot";
            //Start all Files
            File.AppendAllLines(sim_path, new [] {"digraph G {", "start->\"0:\\n" + s.Find(x => x.path == "0").prob + "\";", "start->\"1:\\n" + s.Find(x => x.path == "1").prob + "\";"});
            File.AppendAllLines(num_path, new [] {"digraph G {", "start->\"0:\\n" + n.Find(x => x.path == "0").prob + "\";", "start->\"1:\\n" + n.Find(x => x.path == "1").prob + "\";"});
            File.AppendAllLines(err_path, new [] {"digraph G {", "start->\"0:\\n" + e.Find(x => x.path == "0").prob + "\";", "start->\"1:\\n" + e.Find(x => x.path == "1").prob + "\";"});
            //Write Data Structures into Files
            foreach (PTrace p in s)
            {
                //If error is above threshold then outline node red
                if (e.Find(x => x.path == p.path).prob > 0.001)
                {
                    File.AppendAllLines(sim_path, new [] {"\"" + p.path + ":\\n" + s.Find(x => x.path == p.path).prob + "\" [color=red];"});
                    File.AppendAllLines(num_path, new [] {"\"" + p.path + ":\\n" + n.Find(x => x.path == p.path).prob + "\" [color=red];"});
                    File.AppendAllLines(err_path, new [] {"\"" + p.path + ":\\n" + e.Find(x => x.path == p.path).prob + "\" [color=red];"});
                }
                //If node has a left child add it first
                if (e.Exists(x => x.path == p.path + "0"))
                {
                    File.AppendAllLines(sim_path, new [] {"\"" + p.path + ":\\n" + s.Find(x => x.path == p.path).prob + "\" -> \"" + p.path + "0:\\n" + s.Find(x => x.path == p.path+"0").prob + "\";"});
                    File.AppendAllLines(num_path, new [] {"\"" + p.path + ":\\n" + n.Find(x => x.path == p.path).prob + "\" -> \"" + p.path + "0:\\n" + n.Find(x => x.path == p.path+"0").prob + "\";"});
                    File.AppendAllLines(err_path, new [] {"\"" + p.path + ":\\n" + e.Find(x => x.path == p.path).prob + "\" -> \"" + p.path + "0:\\n" + e.Find(x => x.path == p.path+"0").prob + "\";"});
                }
                //If node has a right child add it next
                if (e.Exists(x => x.path == p.path + "1"))
                {
                    File.AppendAllLines(sim_path, new [] {"\"" + p.path + ":\\n" + s.Find(x => x.path == p.path).prob + "\" -> \"" + p.path + "1:\\n" + s.Find(x => x.path == p.path+"1").prob + "\";"});
                    File.AppendAllLines(num_path, new [] {"\"" + p.path + ":\\n" + n.Find(x => x.path == p.path).prob + "\" -> \"" + p.path + "1:\\n" + n.Find(x => x.path == p.path+"1").prob + "\";"});
                    File.AppendAllLines(err_path, new [] {"\"" + p.path + ":\\n" + e.Find(x => x.path == p.path).prob + "\" -> \"" + p.path + "1:\\n" + e.Find(x => x.path == p.path+"1").prob + "\";"});
                }

            }

            //End all Files
            File.AppendAllLines(sim_path, new [] {"}"});
            File.AppendAllLines(num_path, new [] {"}"});
            File.AppendAllLines(err_path, new [] {"}"});
        } 

        /**** EXTENDED PROBABILITY FUNCTIONS ****/

        //Function which randomly creates a number of components based on the min and max reliability given
        static void GenerateComponents(double min, double max)
        {
            List<double> temp = new List<double>();
            Random R = new Random();
            for (int i = 0; i < number; i++)
                temp.Add(R.NextDouble() * (max-min) + min);
            component_reliabilities = temp.ToArray();
        }
        //Generates Correlations based on component values
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
        //Generates Sigma values based on Correlations
        static void GenerateSigma()
        {
            sigma.Clear();
            foreach (double d in component_reliabilities)
            {
                sigma.Add(Math.Sqrt(d * (1 - d)));
            };
        }

        //Function which calculates bridge reliability based on n number of components and cutsets CS
        static double FindReliability()
        {
            Generate_Bm();
            //Enumerate All Successful States
            List<int[]> Paths = RecursiveEnumerate();
            double R = 0;
            //Loop through and see which Leaf Nodes are success
            foreach (int[] i in Paths)
            {
                    /*Console.Write("Pr(");
                    foreach (int ii in i)
                        Console.Write(ii);
                    Console.WriteLine(") \t= " + TotalProbability(i).ToString("#0.0000000000"));*/
                    R += TotalProbability(i);
                    //Console.WriteLine("Reliability \t= " + R.ToString("#0.0000000000") + "\n");
            }
            Console.WriteLine("Reliability \t= " + R.ToString("#0.0000000000"));
            return R;
        }
        //Recursive Dynamic Approach to Enumerating all possible states, defaults to only adding Full Path and not partial
        static List<int[]> RecursiveEnumerate()
        {
            int[] counters = new int[number];
            for (int i = 0; i < number; i++)
                counters[i] = 0;
            List<int[]> Paths = new List<int[]>();

            Paths = RecursiveEnumerate(0, counters, Paths);
            return Paths;          
        }
        //Second Form of Recursive Function
        static List<int[]> RecursiveEnumerate(int level, int[] counters, List<int[]> Paths)
        {
            List<int> temp_state = new List<int>();
            bool pathflag;
            for (int i = 1; i >= 0; i--)
            {
                temp_state.Clear();
                pathflag = false;
                counters[level] = i;
                for (int a = 0; a <= level; a++)
                    temp_state.Add(counters[a]);
                if ( (level == counters.Length - 1 && !ContainsCutset(temp_state.ToArray())) || ContainsPathset(temp_state.ToArray()))
                {
                    Paths.Add(temp_state.ToArray());
                    pathflag = true;
                }
                if (level != counters.Length - 1 && !ContainsCutset(temp_state.ToArray()) && !pathflag)
                    Paths = RecursiveEnumerate(level + 1, counters, Paths);
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
        //Check to see if the integer array is a pathset, returns true if the query is a pathset, returns false otherwise
        static bool ContainsPathset(int[] query)
        {
            bool safe;
            foreach (int[] path in PS)            //Loop through each cutset and compare to query
            {
                safe = true;
                if (query.Count() >= path.Max())     //See if check is applicable
                {
                    foreach (int i in path)
                    {
                        if (query[i - 1] != 1)
                            safe = false;
                    }
                    if (safe)
                        return true;
                }
                    
            }
            return false;
        }
 
        /**** CORE PROBABILITY FUNCTIONS ****/
        
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

            if (R > 1 || R < 0)
            {
                if(!TotalPErrors.Exists(x => x.path == String.Join("", T1)))
                    TotalPErrors.Add(new PTrace(String.Join("", T1), R));
                Console.WriteLine("Total Probability Error: " + String.Join("", T1) + " = " + R);
            }

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
                    R += (1 - component_reliabilities[i]) * Bm[B.Length][i + 1];
                    if (print)
                        Console.Write("(1 - u[" + (i+1) + "])b(" + (B.Length) + ", " + (i + 1) + ")");
                }
                else
                {
                    R += (0 - component_reliabilities[i]) * Bm[B.Length][i + 1];
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
            {
                if (R > 1 || R < 0)
                {
                    if(!ConditionalPErrors.Exists(x => x.path == String.Join("", B)))
                        ConditionalPErrors.Add(new PTrace(String.Join("", B), R));
                    Console.WriteLine("Conditional Probability Error: " + String.Join("", B) + " = " + R);
                }
                return R;
            }
            else
            {
                if ((1-R) > 1 || (1-R) < 0)
                {
                    if(!ConditionalPErrors.Exists(x => x.path == B.ToString()))
                        ConditionalPErrors.Add(new PTrace(B.ToString(), 1-R));
                    Console.WriteLine("Conditional Probability Error: " + B.ToString() + " = " + (1-R));
                }
                return 1 - R;
            }
        }

        //The b function which returns the appropriate b vector based on level
        static double[] b(int i, bool print = false)
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
            //return R[j-1,0];
            return R.AsColumnMajorArray();
            
        }

        /**** MISC FUNCTIONS ****/

        //Function that should reproduce the results of the RAMS2014 Paper (Make sure correlations matrix and component reliabilities are changed)
        static void RAMS2014()
        {
            //static double[] component_reliabilities = { 0.99, 0.995, 0.98, 0.995};        //RAMS2014
            /*double[] correlations = {1.0000, 0.0094, -0.0023, -0.0020,              //RAMS2014
                                     0.0094, 1.0000,  0.0011,  0.0010,
                                    -0.0023, 0.0011,  1.0000,  0.0010,
                                    -0.0020, 0.0010,  0.0010,  1.0000}; */
            List<int[]> All_States = RecursiveEnumerate();
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

        //Function which checks correlation_matrix to see if the values are within the minimum and maximum bounds set by the equations below
        static List<CorrelationPairs> CorrelationCheck()
        {
            List<CorrelationPairs> bounds = CorrelationBounds();
            List<CorrelationPairs> problems = new List<CorrelationPairs>();

            for (int i = 0; i < component_reliabilities.Count() - 1; i++)
            {
                for (int j = i + 1; j < component_reliabilities.Count(); j++)
                {
                    if ((correlations_matrix[i,j] < bounds.Find(z => z.x == i && z.y == j).min) || (correlations_matrix[i,j] > bounds.Find(z => z.x == i && z.y == j).max))
                        problems.Add(bounds.Find(z => z.x == i && z.y == j));
                }
            }
            return problems;
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
                    */temp.Add(new CorrelationPairs(i, j, min, max, bound));
                }
            }
            return temp;
        }
        //Function which examines multiple runs of Convergence Test
        static void PlotConvergence()
        {
            //success_percentages[0] correlates to the first set of data from m = 0 -> 1.0
            List<double[]> success_percentages = new List<double[]>();
            for (number = 10; number <= 10; number++)
                success_percentages.Add(ConvergenceTest(0.1,1000));
            Console.WriteLine("All Percentages Added");


        }
        //Examines effect of correlation multiplier on convergence of bounds 
        static double[] ConvergenceTest(double step_size, int runs)
        {
            List<double> percentages = new List<double>(); //List to hold percentages of all steps from 0 -> 1.0
            double result;      //resulting percentage after x amount of runs
            double sum;           //Total Reliability of the system, if this value exceeds 1 or is lower than 0, we know there was an error
            bool neg_flag;      //Flag that determines whether or not a set of probabilities had a negative value contained within it
            int success_counter;//Counts the number of successful operations
            
            //Enumerate Final States
            List<int[]> Final_States = RecursiveEnumerate();

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
        //Structure used in multiple functions
        struct CorrelationPairs
        {
            public int x;
            public int y;
            public double min;
            public double max;
            public double value;

            public CorrelationPairs(int xx, int yy, double lower, double upper, double bound = 0)
            {
                x = xx;
                y = yy;
                min = lower;
                max = upper;
                value = bound;
            }
        };

        public class PTrace
        {
            public string path;
            public int count;
            public double prob;

            public PTrace(string p, double pr = 0)
            {
                path = p;
                count = 1;
                prob = pr;
            }
        };

        /**** LEGACY FUNCTIONS ****/

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

            All_States = RecursiveEnumerate();
            
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
        //Prints out Total Probability for all 32 Combinations
        static void PrintTotal32()
        {
            //Print out Component Reliabilites and Correlations Matrix
            Console.Write("Component Reliabilities: {");
            component_reliabilities.ToList().ForEach(i => Console.Write(i.ToString() + ", "));
            Console.WriteLine("}\n");
            Console.WriteLine("Correlations Matrix:");
            Console.WriteLine(correlations_matrix.ToString());
            List<int[]> FinalStates = RecursiveEnumerate();
            foreach (int[] i in FinalStates)
            {
                Console.Write("Pr(");
                foreach (int j in i)
                    Console.Write(j);
                Console.WriteLine(")\t=\t" + TotalProbability(i).ToString("#0.############"));
            }
         }

    }
}
