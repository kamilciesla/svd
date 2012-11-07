using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using System.Reflection;
using SVD;

namespace svd
{

    public class Program
    {
        public static void MatrixProduct(float[,] result, float[,] Left, float[,] Right, int nRows, int nCommon, int nCols)
        {
            for (int i = 0; i < nRows; ++i)
            {
                for (int j = 0; j < nCols; ++j)
                {
                    result[i, j] = (float)(0.0);
                    for (int k = 0; k < nCommon; ++k)
                    {
                        result[i, j] += Left[i, k] * Right[k, j];
                    }
                }
            }
        }

        public static void PrintMatrix(string name, float[,] Matrix, int nRows, int nCols, bool revert = false)
        {

            Console.WriteLine(name);
            if (!revert)
            {
                for (int i = 0; i < nRows; ++i)
                {
                    for (int j = 0; j < nCols; ++j)
                    {
                        Console.Write(" {0:###0.00} ", Matrix[i, j]);
                    }
                    Console.WriteLine();
                }
            }
            else
            {
                for (int i = 0; i < nCols; ++i)
                {
                    for (int j = 0; j < nRows; ++j)
                    {
                        Console.Write(" {0:###0.00} ", Matrix[j, i]);
                    }
                    Console.WriteLine();
                }

            }
            Console.WriteLine();
        }


        //dane wejściowe
        static float[,] matrix = new float[,]
        {
            {1f,0f,0f,1f,0f,0f,0f,0f,0f},
            {1f,0f,1f,0f,0f,0f,0f,0f,0f},
            {1f,1f,0f,0f,0f,0f,0f,0f,0f},
            {0f,1f,1f,0f,1f,0f,0f,0f,0f},
            {0f,1f,1f,2f,0f,0f,0f,0f,0f},
            {0f,1f,0f,0f,1f,0f,0f,0f,0f},
            {0f,1f,0f,0f,1f,0f,0f,0f,0f},
            {0f,0f,1f,1f,0f,0f,0f,0f,0f},
            {0f,1f,0f,0f,0f,0f,0f,0f,1f},
            {0f,0f,0f,0f,0f,1f,1f,1f,0f},
            {0f,0f,0f,0f,0f,0f,1f,1f,1f},
            {0f,0f,0f,0f,0f,0f,0f,1f,1f}
        };


        public static float skalarny(float[] wektor1, float[] wektor2)
        {
            float skalar = 0;

            float[] tmp = new float[wektor1.Length];
            for (int i = 0; i < wektor1.Length; i++)
            {
                tmp[i] = wektor1[i] * wektor2[i];
                skalar += tmp[i];
            }
            return skalar;
        }
        public static double dlugoscVectora(float[] wektor)
        {
            float dlugosc = 0;

            for (int i = 0; i < wektor.Length; i++)
            {
                dlugosc += wektor[i] * wektor[i];
            }

            return Math.Sqrt(dlugosc);
        }

        public static float[] getVector(float[,] inputMatrix, int nCols, int index)
        {
            float[] vector = new float[nCols];

            for (int i = 0; i < nCols; i++)
            {
                vector[i] = inputMatrix[index, i];
            }

            return vector;
        }

        static double dotProduct = 0;
        private static void execute()
        {
            float[,] res = new float[8, 8];

            float[] vector1 = new float[9];
            float[] vector2 = new float[9];

            int index = 0;
            for (int x = 0; x < 9; x++)
            {
                vector1 = getVector(UxSxVmatrix, 9, index);
                for (int y = 1 + index; y < 12; y++)
                {
                    vector2 = getVector(UxSxVmatrix, 9, y);
                    float skalar = skalarny(vector1, vector2);
                    double dl_v1 = dlugoscVectora(vector1);
                    double dl_v2 = dlugoscVectora(vector2);

                    dotProduct = (double)skalar/(double)(dl_v1*dl_v2);
                    try
                    {
                        res[y-1, x] = (float)dotProduct;
                    }
                    catch { }
                }
                index++;
            }

            PrintMatrix("after svd", res, 8, 8);
        }

        static int nCols = 9;
        static int nRows = 12;

        static float[,] U;
        static float[,] U2;
        
        static float[,] V;
        static float[,] V2;

        static float[,] S;
        static float[,] S2;

        static float[,] UxSmatrix;
        static float[,] UxSxVmatrix;

        private static void GetS()
        {
            //Read and make matrix with singular values. It is diagonal.
            S = new float[nCols, nRows];
            for (int i = 0; i < nCols; ++i)
            {
                for (int j = 0; j < nRows; ++j)
                {
                    S[i, j] = (float)(0.0);
                }
            }
            int rank = nCols;

            string singularValues = resultBlockName + "-S";
            string line = string.Empty;
            System.IO.StreamReader file = new System.IO.StreamReader(singularValues);
            line = file.ReadLine();
            try
            {
                int cnt = 0;
                while ((line = file.ReadLine()) != null)
                {
                    S[cnt, cnt] = (float)(Convert.ToDouble(line));
                    ++cnt;
                }
            }
            catch (Exception)
            {
                Console.WriteLine("Misformatted file: {0}", singularValues);
                Environment.Exit(1);
            }
            file.Close();
            PrintMatrix("S", S, rank, rank);
        }
        private static void GetS2()
        {
            S2 = new float[2, 2];
            for (int i = 0; i < 2; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    S2[i, j] = S[i, j];
                }
            }

            PrintMatrix("S [2x2]", S2, 2, 2);
        }

        private static void GetU()
        {
            int rank = nCols;

            //Read and make UT matrix, it is transposed U in U*S*VT
            U = new float[nCols, nRows];
            for (int i = 0; i < nCols; ++i)
            {
                for (int j = 0; j < nRows; ++j)
                {
                    U[i, j] = (float)(0.0);
                }
            }
            string fUTfileName = resultBlockName + "-Ut";
            using (FileStream stream = new FileStream(fUTfileName, FileMode.Open))
            {
                using (BinaryReader reader = new BinaryReader(stream))
                {
                    int nRowsRead = System.Net.IPAddress.NetworkToHostOrder(reader.ReadInt32());
                    int nColsRead = System.Net.IPAddress.NetworkToHostOrder(reader.ReadInt32());
                    for (int i = 0; i < nRowsRead; ++i)
                    {
                        for (int j = 0; j < nColsRead; ++j)
                        {
                            int nBuf = System.Net.IPAddress.NetworkToHostOrder(reader.ReadInt32());
                            byte[] b = BitConverter.GetBytes(nBuf);
                            U[i, j] = BitConverter.ToSingle(b, 0);
                        }
                    }
                    reader.Close();
                }
                stream.Close();
            }
            //end reading


            PrintMatrix("\nU", U, rank, nRows, true);
        }
        private static void GetU2()
        {
            U2 = new float[12, 2];
            for (int i = 0; i < 2; ++i)
            {
                for (int j = 0; j < 12; ++j)
                {
                    U2[j, i] = U[i, j];
                }
            }

            PrintMatrix("U [12x2]", U2, 12, 2);
        }




        private static void GetV()
        {
            int rank = nCols;
            //Read and make VT matrix
            V = new float[nCols, nCols];
            for (int i = 0; i < nCols; ++i)
            {
                for (int j = 0; j < nCols; ++j)
                {
                    V[i, j] = (float)(0.0);
                }
            }
            string fVTfileName = resultBlockName + "-Vt";
            using (FileStream stream = new FileStream(fVTfileName, FileMode.Open))
            {
                using (BinaryReader reader = new BinaryReader(stream))
                {
                    int nRowsRead = System.Net.IPAddress.NetworkToHostOrder(reader.ReadInt32());
                    int nColsRead = System.Net.IPAddress.NetworkToHostOrder(reader.ReadInt32());
                    for (int i = 0; i < nRowsRead; ++i)
                    {
                        for (int j = 0; j < nColsRead; ++j)
                        {
                            int nBuf = System.Net.IPAddress.NetworkToHostOrder(reader.ReadInt32());
                            byte[] b = BitConverter.GetBytes(nBuf);
                            V[i, j] = BitConverter.ToSingle(b, 0);
                        }
                    }
                    reader.Close();
                }
                stream.Close();
            }
            //end reading


            PrintMatrix("\nV", V, rank, rank);
        }
        private static void GetV2()
        {
            V2 = new float[2, 9];
            for (int i = 0; i < 2; ++i)
            {
                for (int j = 0; j < 9; ++j)
                {
                    V2[i, j] = V[j, i];
                }
            }

            PrintMatrix("V [2x9]", V2, 2, 9);
        }

        private static void UxS()
        {
            UxSmatrix = new float[12,2];

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 12; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        UxSmatrix[j, i] += U[k, j] * S[k, i];
                    }
                }
            }

            PrintMatrix("UxS", UxSmatrix, 12, 2);
        }

        private static void UxSxV()
        {
            UxSxVmatrix = new float[12, 9];

            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 12; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        UxSxVmatrix[j, i] += UxSmatrix[j, k] * V[k, i];
                    }
                }
            }

            PrintMatrix("UxS xV", UxSxVmatrix, 12, 9);
        }

        static string dataFileName;
        static string resultBlockName;

        static void Main(string[] args)
        {
            const double kappa = 1e-6;
            string path = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
            dataFileName = Path.Combine(path, "matrix");
            resultBlockName = Path.Combine(path, "result");




            #region Write input matrix to file
            int nNonZeros = 0;
            int[] nNonZerosInCols = new int[nCols];
            for (int i = 0; i < nCols; i++)
            {
                nNonZerosInCols[i] = 0;
            }
            for (int j = 0; j < nCols; ++j)
            {
                for (int i = 0; i < nRows; ++i)
                {
                    if (Math.Abs(matrix[i, j]) > kappa)
                    {
                        ++nNonZerosInCols[j];
                        ++nNonZeros;
                    }
                }
            }


            using (FileStream stream = new FileStream(dataFileName, FileMode.Create))
            {
                using (BinaryWriter writer = new BinaryWriter(stream))
                {
                    writer.Write(System.Net.IPAddress.HostToNetworkOrder(nRows));
                    writer.Write(System.Net.IPAddress.HostToNetworkOrder(nCols));
                    writer.Write(System.Net.IPAddress.HostToNetworkOrder(nNonZeros));
                    for (int k = 0; k < nCols; ++k)
                    {
                        writer.Write(System.Net.IPAddress.HostToNetworkOrder(nNonZerosInCols[k]));
                        for (int j = 0; j < nRows; ++j)
                        {
                            if (Math.Abs(matrix[j, k]) > kappa)
                            {
                                writer.Write(System.Net.IPAddress.HostToNetworkOrder(j));
                                byte[] b = BitConverter.GetBytes(matrix[j, k]);
                                int x = BitConverter.ToInt32(b, 0);
                                writer.Write(System.Net.IPAddress.HostToNetworkOrder(x));
                            }
                        }
                    }
                    writer.Flush();
                    writer.Close();
                }
                stream.Close();
            }
            #endregion

            SVD.SVD.ProcessData(dataFileName, resultBlockName);

            PrintMatrix("\nInput matrix", matrix, nRows, nCols);

            GetU();
            GetV();
            GetS();

            GetU2();
            GetV2();
            GetS2();

            UxS();

            UxSxV();

            execute();

            Console.WriteLine("Press any key to stop");
            do
            {
                while (!Console.KeyAvailable)
                {
                    // Do something
                }
            } while (Console.ReadKey(true).Key != ConsoleKey.Escape);

        }
    }
}
