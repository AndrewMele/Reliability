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

            CoverageTest2(10000);
            Console.ReadLine();
            //Build environment
            Build("3o5");
            
            /*/print info
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
            Console.WriteLine("");*/
            
            //Execute run
            //FindReliability();
            //List<CorrelationPairs> issues = CorrelationCheck();
            //Console.WriteLine(correlations_matrix.ToString());

            FindReliability();
            double S = Simulation(1000000);
            Console.WriteLine("Simulation Estimate with 1,000 runs: " + S.ToString("#0.0000000000"));

            //Hold Results on screen
            endTime = DateTime.Now;
            elpasedMilliseconds = ((TimeSpan)(endTime - startTime)).TotalMilliseconds;
            Console.WriteLine("Task took " + elpasedMilliseconds.ToString("#0.00") + " ms.");           
        }
        //Functions that need creating... Setup function: switch case to setup environment for each test case
        //                                Run function: switch case to run appropriate example for each test case
        
        //First Coverage Test Function which tests a 2o3 enviornment with random reliabilities and correlations within the theoretical bounds
        static void CoverageTest1(int runs)
        {
            Random r = new Random();        //Random Number Generator to generate reliabilities
            int success = 0;                //Initialize success counter to zero
            number = 3;                     //Number of components to generate
            CS.Clear();                     //Clear any previous cutsets stored
            PS.Clear();                     //Clear any previous pathsets stored
            //Add cuts and paths for a 2o3 component system
            CS.Add(new int[] {1,2});
            CS.Add(new int[] {1,3});
            CS.Add(new int[] {2,3});
            PS.Add(new int[] {1,2});
            PS.Add(new int[] {1,3});
            PS.Add(new int[] {2,3});

            for (int i = 0; i < runs; i++)  //Loop for the appropriate number of test runs 
            {
                Bm.Clear();                     //Clear B matrix and reinitialize the vector with 2 placeholders
                Bm.Add(new double[] {0});
                Bm.Add(new double[] {0});
                GenerateComponents(0.01, 0.99); //Generate 3 random reliabilities
                GenerateCorrelations();         //Generate correlations within the theoretical bounds
                GenerateSigma();                //Generate Sigma values for use in B matrix

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
        static void CoverageTest2(int runs)
        {
            Random r = new Random();        //Random Number Generator to choose random correlations from the data set
            int switch_index = 0;           //Variable used in correlation assurance, 0 = c01, 1 = c02, 2 = c12
            int a;                          //Temp variable to hold random number between 1-162 (used in correlation generation)
            double b1,b2;                   //Temp variables to hold bounds for random number generation when generating correlations
            double c01 = 0,c02 = 0,c12 = 0; //Temp variables to hold potential correlations while we perform tests 
                                            //c01 = correlation between components 1 and 2  c02 = correlation between components 1 and 3  c12 = correlation between components 2 and 3
            double[] correlations;          //Temp variable used in correlation matrix construction
            int success = 0;                //Initialize success counter to zero
            number = 3;                     //Number of components to generate
            CS.Clear();                     //Clear any previous cutsets stored
            PS.Clear();                     //Clear any previous pathsets stored
            //Add cuts and paths for a 2o3 component system
            CS.Add(new int[] {1,2});
            CS.Add(new int[] {1,3});
            CS.Add(new int[] {2,3});
            PS.Add(new int[] {1,2});
            PS.Add(new int[] {1,3});
            PS.Add(new int[] {2,3});
            //Create Beta Distribution with alpha = 1000 and beta = 2, add 10,000 samples to double array and compute average
            Beta beta = new Beta(1000,2);
            double[] samples = new double[3*runs];
            beta.Samples(samples);
            double average = samples.Average(); 
            //Make sure average of beta distribution is within the bounds
            if (average < 0.997 || average > 0.999)
            {
                Console.WriteLine("Beta Distribution Error: Average sample = " + average.ToString("#.######"));
                return;
            }

            //Correlation data set
            double[] data_set = new double[] {0.00617284, 0.0123457, 0.0185185, 0.0246914, 0.0308642, 0.037037, 0.0432099, 0.0493827, 0.0555556, 0.0617284, 0.0679012, 0.0740741, 0.0802469, 0.0864198, 0.0925926, 0.0987654, 0.104938, 0.111111, 0.117284, 0.123457, 0.12963, 0.135802, 0.141975, 0.148148, 0.154321, 0.160494, 0.166667, 0.17284, 0.179012, 0.185185, 0.191358, 0.197531, 0.203704, 0.209877, 0.216049, 0.222222, 0.228395, 0.234568, 0.240741, 0.246914, 0.253086, 0.259259, 0.265432, 0.271605, 0.277778, 0.283951, 0.290123, 0.296296, 0.302469, 0.308642, 0.314815, 0.320988, 0.32716, 0.333333, 0.339506, 0.345679, 0.351852, 0.358025, 0.364198, 0.37037, 0.376543, 0.382716, 0.388889, 0.395062, 0.401235, 0.407407, 0.41358, 0.419753, 0.425926, 0.432099, 0.438272, 0.444444, 0.450617, 0.45679, 0.462963, 0.469136, 0.475309, 0.481481, 0.487654, 0.493827, 0.5, 0.506173, 0.512346, 0.518519, 0.524691, 0.530864, 0.537037, 0.54321, 0.549383, 0.555556, 0.561728, 0.567901, 0.574074, 0.580247, 0.58642, 0.592593, 0.598765, 0.604938, 0.611111, 0.617284, 0.623457, 0.62963, 0.635802, 0.641975, 0.648148, 0.654321, 0.660494, 0.666667, 0.67284, 0.679012, 0.685185, 0.691358, 0.697531, 0.703704, 0.709877, 0.716049, 0.722222, 0.728395, 0.734568, 0.740741, 0.746914, 0.753086, 0.759259, 0.765432, 0.771605, 0.777778, 0.783951, 0.790123, 0.796296, 0.802469, 0.808642, 0.814815, 0.820988, 0.82716, 0.833333, 0.839506, 0.845679, 0.851852, 0.858025, 0.864198, 0.87037, 0.876543, 0.882716, 0.888889, 0.895062, 0.901235, 0.907407, 0.91358, 0.919753, 0.925926, 0.932099, 0.938272, 0.944444, 0.950617, 0.95679, 0.962963, 0.969136, 0.975309, 0.981481, 0.987654, 0.993827, 1};

            for (int i = 0; i < runs; i++)  //Loop for the appropriate number of test runs 
            {
                Bm.Clear();                     //Clear B matrix and reinitialize the vector with 2 placeholders
                Bm.Add(new double[] {0});
                Bm.Add(new double[] {0});
                //Add 3 samples as component reliabilities then generate correlation bounds
                component_reliabilities = new double[] {samples[i], samples[i+1], samples[i+2]};
                List<CorrelationPairs> bounds = CorrelationBounds();
                //Continue generating correlations using data set as reference until all three correlations are within theoretical bounds
                while(true)
                {
                    if (switch_index == 3)  //If all 3 correlations are within the bounds reset switch index and break loop
                    {
                        switch_index = 0;
                        break;
                    }

                    a = r.Next(1,162);      //Generate random number between 1-162
                    if (a != 161)           //If index does not overflow set b1 and b2 using increment
                    {
                        b1 = data_set[a];
                        b2 = data_set[a+1];
                    }
                    else                    //Else set b1 and b2 using decrement
                    {
                        b1 = data_set[a-1];
                        b2 = data_set[a];
                    }
                    //Switch case statement to determine which correlation we are currently checking, 0 = 12  1 = 13  2 = 23
                    switch (switch_index)
                    {
                        case 0:
                            c01 = r.NextDouble() * (b2-b1) + b1;    //Generate random double between b1 and b2
                            //If c01 is within the bounds increment index otherwise repeat
                            if ((c01 > bounds.Find(z => z.x == 0 && z.y == 1).min) && (c01 < bounds.Find(z => z.x == 0 && z.y == 1).max))
                                switch_index++;
                            break;
                        case 1:
                            c02 = r.NextDouble() * (b2-b1) + b1;    //Generate random double between b1 and b2
                            //If c02 is within the bounds increment index otherwise repeat
                            if ((c02 > bounds.Find(z => z.x == 0 && z.y == 2).min) && (c02 < bounds.Find(z => z.x == 0 && z.y == 2).max))
                                switch_index++;
                            break;
                        case 2:
                            c12 = r.NextDouble() * (b2-b1) + b1;    //Generate random double between b1 and b2
                            //If c12 is within the bounds increment index otherwise repeat
                            if ((c12 > bounds.Find(z => z.x == 1 && z.y == 2).min) && (c12 < bounds.Find(z => z.x == 1 && z.y == 2).max))
                                switch_index++;
                            break;
                        default:
                            Console.WriteLine("Error switch index = " + switch_index);
                            return;
                    }
                }
                
                //Enter correlations into a double array then generate matrix
                correlations = new double[] { 1.00, c01, c02,
                                              c01, 1.00, c12,
                                              c02, c12, 1.00};
                correlations_matrix = Matrix<double>.Build.Dense(3,3,correlations);
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
            Console.WriteLine("Coverage Test 2:\n" + success + " successful enumerations after " + runs + " test runs");
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
                    correlations = new double[] { 1.00, -0.28,
                                                -0.28,  1.00};
                    correlations_matrix = Matrix<double>.Build.Dense(2,2,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1});
                    CS.Add(new int[] {2});
                    PS.Add(new int[] {1,2});
                    break;
                case "2p":                      //2 component parallel
                    number = 2;
                    component_reliabilities = new double[] {0.8, 0.75};
                    correlations = new double[] { 1.00, -0.28,
                                                -0.28,  1.00};
                    correlations_matrix = Matrix<double>.Build.Dense(2,2,correlations);
                    GenerateSigma();
                    CS.Add(new int[] {1,2});
                    PS.Add(new int[] {1});
                    PS.Add(new int[] {2});
                    break;
                case "3s":                      //3 component series
                    number = 3;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7};
                    correlations = new double[] { 1.00, -0.28, 0.1,
                                                -0.28,  1.00, 0.2,
                                                0.1,   0.2,  1.00};
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
                    correlations = new double[] { 1.00, -0.28, 0.1,
                                                -0.28,  1.00, 0.2,
                                                0.1,   0.2,  1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15,
                                                 -0.28,  1.00,   0.2, 0.30,
                                                   0.1,   0.2,  1.00, 0.32,
                                                 -0.15,  0.30,  0.32, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15,
                                                 -0.28,  1.00,   0.2, 0.30,
                                                   0.1,   0.2,  1.00, 0.32,
                                                 -0.15,  0.30,  0.32, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23, 0.22,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,-0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17, 0.07,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 0.25,
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16, 0.10,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,-0.07,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00, 0.35,
                                                  0.22, -0.16,  0.07, 0.25, 0.10,-0.07, 0.35, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15, 0.05,-0.23, 0.22,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,-0.20, 0.16,-0.16,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,-0.30, 0.17, 0.07,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22, 0.20,-0.30, 0.25,
                                                  0.15, -0.15,  0.50,-0.22, 1.00, 0.25,-0.16, 0.10,
                                                  0.05, -0.20, -0.30, 0.20, 0.25, 1.00, 0.27,-0.07,
                                                 -0.23,  0.16,  0.17,-0.30,-0.16, 0.27, 1.00, 0.35,
                                                  0.22, -0.16,  0.07, 0.25, 0.10,-0.07, 0.35, 1.00};
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
                    correlations = new double[] { 1.00, -0.28,   0.1,-0.15, 0.15,
                                                 -0.28,  1.00,   0.2, 0.30,-0.15,
                                                   0.1,   0.2,  1.00, 0.32, 0.50,
                                                 -0.15,  0.30,  0.32, 1.00,-0.22,
                                                  0.15, -0.15,  0.50,-0.22, 1.00};
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
                case "2o3":
                    number = 3;
                    component_reliabilities = new double[] {0.8, 0.75, 0.7};
                    correlations = new double[] { 1.00, -0.28, 0.1,
                                                -0.28,  1.00, 0.2,
                                                0.1,   0.2,  1.00};
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
                    //correlations_matrix = Matrix<double>.Build.Dense(7,7,correlations);
                    GenerateCorrelations(0.2);
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
            //Generate_Bm();
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
            foreach (PTrace p in sim)
            {
                p.prob = (double) p.count / (double) runs;
                num.Add(new PTrace(p.path, TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c)))));
                err.Add(new PTrace(p.path, Math.Abs(p.prob - TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c))))));
                
                //if (Math.Abs(p.prob - TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c)))) > 0.001)
                   // totals.Add(new PTrace(p.path, Math.Abs(p.prob - TotalProbability(Array.ConvertAll(p.path.ToArray(), c => (int) Char.GetNumericValue(c))))));
            }
            CreateDOTGraphs(sim,num,err);
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
            //Console.WriteLine("Reliability \t= " + R.ToString("#0.0000000000"));
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
