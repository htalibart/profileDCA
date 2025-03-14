#ifndef DP_MAT_APURVA_1
#define DP_MAT_APURVA_1

using namespace std;

/**
* DP_Mat_Apurva.
* Dynamic Programming used by the A_purva model for Contact Map Overlap Maximization Problem (CMO).
*/
class dp_mat_apurva : public dp_mat
{
    protected :

    public :

        /**
        * Constructors
        */
        dp_mat_apurva(graph_apurva & g, float * insert_open_row_, float * insert_open_col_, float * insert_extend_row_, float * insert_extend_col_): dp_mat()
            {
        		//gap_open = go;
        		//gap_extend = ge;

		insert_open_row = insert_open_row_;
		insert_open_col = insert_open_col_;
		insert_extend_row = insert_extend_row_;
		insert_extend_col = insert_extend_col_;



                nb_col = g.get_nb_col();
                nb_row = g.get_nb_row();

                dp_arc = new double* [2000];
                arc_from = new int* [2000];
                for(int i(0); i<2000;++i)
                {
                    dp_arc[i] = new double [2000];
                    fill_n( dp_arc[i], 2000, static_cast<double>( 0. ) );

                    arc_from[i] = new int [2000];
                    fill_n( arc_from[i], 2000, static_cast<int>( 0 ) );
                }
                fill_n( arc_from[0], 2000, static_cast<int>( -1 ) );

                dp_arc_out_col = new int** [nb_col];
                for (int i(0); i != nb_col; ++i)
                {
                    dp_arc_out_col[i] = new int* [nb_row];
                    for (int j(0); j != nb_row; ++j)
                    {
                        dp_arc_out_col[i][j] = new int [g.get_nb_next_col(i,j)];
                        fill_n(dp_arc_out_col[i][j], g.get_nb_next_col(i,j), static_cast<int>( -1 ) );
                    }
                }

                dp_M = new double* [nb_col+1];
                dp_GA = new double* [nb_col+1];
                dp_GB = new double* [nb_col+1];
                dp_score = new double* [nb_col+1];
                dp_M_from = new int* [nb_col+1];
                dp_GA_from = new int* [nb_col+1];
                dp_GB_from = new int* [nb_col+1];
                for (int i(0); i != nb_col+1; ++i)
                {
                    dp_M[i] = new double [nb_row+1];
                    //fill_n(dp[i], nb_row+1, static_cast<double>(0.));

                    dp_GA[i] = new double [nb_row+1];
                    //fill_n(dp_v[i], nb_row+1, static_cast<double>(0.));

                    dp_GB[i] = new double [nb_row+1];
                    //fill_n(dp_h[i], nb_row+1, static_cast<double>(0.));

                    dp_score[i] = new double [nb_row+1];
                    fill_n(dp_score[i], nb_row+1, static_cast<double>(0.));

                    dp_M_from[i] = new int [nb_row+1];
                    //fill_n(dp_from[i], nb_row+1, static_cast<int>(0));

                    dp_GA_from[i] = new int [nb_row+1];
                    //fill_n(dp_v_from[i], nb_row+1, static_cast<int>(0));

                    dp_GB_from[i] = new int [nb_row+1];
                    //fill_n(dp_h_from[i], nb_row+1, static_cast<int>(0));
                }

                //Initialization of DP matrices
                //diagonal: dp (symbols aligned)
                //horizontal: dp_h (gap in first sequence), Q
                //vertical: dp_v (gap in second sequence), P

                dp_M[0][0] = 0.0;
                dp_GB[0][0] = 0.0;
                dp_GA[0][0] = 0.0;

                for(int i=1; i<=nb_row; ++i)
                {
			dp_M[0][i]=INFINITY;
			dp_GA[0][i]=INFINITY;
			dp_GB[0][i]=get_insertion_open_after_col(-1)+(i-1)*get_insertion_extend_after_col(-1);
                }


                for(int k=1; k<=nb_col; ++k)
                {
			dp_M[k][0]=INFINITY;
			dp_GB[k][0]=INFINITY;
			dp_GA[k][0]=get_insertion_open_after_row(-1)+(k-1)*get_insertion_extend_after_row(-1);
                }

            };

        /**
        * Destructor
        */
        ~dp_mat_apurva()
        {
            for(int i(0); i<2000;++i)
            {
                delete [] dp_arc[i];
                delete [] arc_from[i];
            }
            delete [] dp_arc;
            delete [] arc_from;

            for (int i(0); i != nb_col; ++i)
            {
                for (int j(0); j != nb_row; ++j)
                {
                    delete [] dp_arc_out_col[i][j];
                }
                delete [] dp_arc_out_col[i];
            }
            delete [] dp_arc_out_col;

            for (int i(0); i != nb_col+1; ++i)
            {
                delete [] dp_M[i];
                delete [] dp_score[i];
                delete [] dp_M_from[i];
                delete [] dp_GA[i];
                delete [] dp_GA_from[i];
                delete [] dp_GB[i];
                delete [] dp_GB_from[i];
            }
            delete [] dp_M;
            delete [] dp_score;
            delete [] dp_M_from;
            delete [] dp_GA;
            delete [] dp_GA_from;
            delete [] dp_GB;
            delete [] dp_GB_from;
        };

        /**
        * Methods
        */

        /**
        * Fill the dp_score table, ie compute for each node, the value of outgoing edges set.
        */
        void fill(graph_apurva & g, lambda_mat_apurva & lb_mat, int * lo, int * up);
        /**
        * Find the node set having the best outgoing edges values
        */
        double solve(graph_apurva & g, int * sol, lambda_mat_apurva &, int * lo, int * up);
        double solve_w_gapcosts(graph_apurva & g, int * sol, int * sol_insert_before, lambda_mat_apurva & lb_mat, int * lo, int * up);

        /**
        * For a given node N_{col1,row1} and a given next column (define by a next_col indice ind), return
        * the chosen row2, corresponding to the activated edge E_{col1,row1,col2,row2}.
        */
        inline int get_arc_out(int col1, int row1, int ind){ return dp_arc_out_col[col1][row1][ind]; }
        /**
        * For a given node, find the best value of the best outgoing edges set.
        */
        double value_arcs_out(int nb_next_col, int nb_next_row, int * next_col, int * next_row, double * coef_lambda_edge);
        /**
        * For a given node, find the best value of the best outgoing edges set, and fill dp_arc_out_col
        */
        void calc_arcs_out(int col1, int row1, graph_apurva & g, lambda_mat_apurva & lb_mat, int * lo, int * up);

        /**
         * Get the gap open and gap extension costs
         */
        //inline double get_gap_open(){ return gap_open; }
        //inline double get_gap_extend(){ return gap_extend; }
	//inline float get_insertion_open_after_col(int col){cout << "after col " << col << " : " << insert_open_col[col+1] << endl; return insert_open_col[col+1];}
	inline float get_insertion_open_after_col(int col)
	{
		if (col<-1 || col>nb_col)
		{
			cout << "insertion_open_after_col" << col << endl;
		}
		return insert_open_col[col+1];
	}
	//inline float get_insertion_open_after_row(int row){cout << "after row " << row << " : " << insert_open_row[row+1] << endl; return insert_open_row[row+1];}
	inline float get_insertion_open_after_row(int row)
	{
		if (row<-1 || row>nb_row)
		{
			cout << "insertion_open_after_row" << row << endl;
		}

		return insert_open_row[row+1];
	}
	inline float get_insertion_extend_after_col(int col)
	{
		if (col<-1 || col>nb_col)
		{
			cout << "insertion_extend_after_col" << col << endl;
		}

		return insert_extend_col[col+1];
	}
	inline float get_insertion_extend_after_row(int row)
	{
		if (row<-1 || row>nb_row)
		{
			cout << "insertion_extend_after_row" << row << endl;
		}

		return insert_extend_row[row+1];
	}



};

#endif
