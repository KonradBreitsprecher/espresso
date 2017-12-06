#include <boost/multi_array.hpp>
#include "grid.hpp"


template <typename T, size_t N>
class LatticeMA //: Communication::InstanceCallback
{
    public:
        //---Typedefs---
        typedef boost::multi_array<T, N> tArray;
        typedef boost::const_multi_array_ref<T, N> tArrayRef;
        typedef typename tArray::index tIndex;
        typedef boost::array<tIndex, N> tIndexArray;
        typedef boost::multi_array_types::index_range range;
        typedef typename tArray::template array_view<N>::type tView;
        
        typedef boost::multi_array< boost::array<T,N> , N> tArrayGradient;

    private:        
        //---Helpers---

        // Negative modulo
        int mod(int val, int mod_val) {
            return ((val % mod_val) + mod_val) % mod_val;
        }
        
        // Flat index to indexArray
        template <typename U>
        tIndexArray getIndexArray( const U& m, int flat_index)
        {
          tIndexArray _index;
          for ( unsigned int dir = 0; dir < N; dir++ ) {
             _index[dir] = (flat_index / m.strides()[dir]) % m.shape()[dir] +  m.index_bases()[dir]; 
          }
          return _index;
        }

        // Index array in local lattice data from position
        bool get_local_index_from_position(double *pos, tIndexArray &local_index)
        {
            for(int dim = 0; dim < N; ++dim) {
                int l = int(pos[dim]/box_l[dim]*m_bins[dim]) - m_local_halo_offset[dim];
                //std::cout << dim << " " << pos[dim] << " " << m_bins[dim] <<  " " <<  l << " " << m_halo_size[dim] << " " << m_local_block_size[dim] << std::endl;
                if (l < -m_halo_size[dim] || l >= m_local_block_size[dim]-m_halo_size[dim]) {
                    return false;
                } else {
                    local_index[dim] = l;
                }
            }
            return true;
        }
        
        // Index array in global lattice data from position
        tIndexArray get_global_index_from_position(double *pos)
        {
            tIndexArray global_index;
            for(int dim = 0; dim < N; ++dim) {
                global_index[dim] = int(pos[dim]/box_l[dim]*m_bins[dim]);
            }
            return global_index; 
        }

        //---Sizing and data import---

        // Calculates and resizes the local array according to the bins per dimension, halo and node position in the node grid.
        // Assigns the remainder to the last node in each dimension if the bins (without halo) cannot be distributed equally on the nodes
        void size_local_data(const int *node_pos, const int *node_grid)
        {
            // Get local array size
            boost::array<size_t, N> extents;
            for(int dim = 0; dim < N; ++dim) {
                //Last node gets remainder, if bins[dim] is not integer multiple of node_grid[dim] 
                //If there are more nodes than bins in this dim: assign single bin to all nodes, limit remainder and offset
                int int_bins_per_node = std::max(1, int((int)m_bins[dim]/node_grid[dim]));
                int remainder = std::max(0,(int)m_bins[dim] - node_grid[dim]*int_bins_per_node);
                int block_size = int_bins_per_node+2*m_halo_size[dim];
            
                if (node_pos[dim] < node_grid[dim]-1) {
                    m_local_block_size[dim] = block_size;
                } else {
                    m_local_block_size[dim] = block_size + remainder ;
                }

                m_local_halo_offset[dim] = std::min((int)m_bins[dim]-1, node_pos[dim]*int_bins_per_node);
                extents[dim] = m_local_block_size[dim];
            }

            //--->Size the local mutiarray = data + halo
            m_local_data.resize(extents);
        }

        // Takes the global data and sets the local data (with halo) 
        void assign_global_to_local_data(tArrayRef data)
        {
            tIndexArray local_index;
            tIndexArray global_index;
            // Run through all elements of local_data...
            for ( int i = 0; i < m_local_data.num_elements(); i++ ) {
                // ...get it's index array.
                local_index = getIndexArray( m_local_data, i );
                for(int dim = 0; dim < N; ++dim) {
                    // Calculate the global index array respecting halos...
                    global_index[dim] = mod( m_local_halo_offset[dim] - m_halo_size[dim] + local_index[dim], m_bins[dim] );
                }
                // ...and assign global to local data.
                m_local_data(local_index) = data(global_index);
            }
        }

        // Rebase local data so that negative indices are halo and the world data starts at 0
        void rebase_local_data()
        {
            boost::array<tIndex,N> bases;  
            for(int dim = 0; dim < N; ++dim) {
                bases[dim] = -m_halo_size[dim];
            }
            m_local_data.reindex(bases);
        }
        

        //---Members---
        std::array<int, N> m_bins;              // bins
        std::array<int, N> m_halo_size;         // halo
        std::array<int, N> m_local_halo_offset; // local data starts here in global data
        std::array<int, N> m_local_block_size;  // local block size (block = data + halo)
        tArray m_local_data;                    // the actual local data block
       
        tArrayGradient m_local_data_gradient;   // precalculated gradient

    public:

        LatticeMA(tArrayRef data, std::array<int, N> halo_size) 
            : m_halo_size(halo_size)
        {
            // Bins from incoming data shape (for global data, without halo)
            for(int dim = 0; dim < N; ++dim) {
                m_bins[dim] = data.shape()[dim];
                assert(m_bins[dim] > 0);
                assert(m_halo_size[dim] >= 0);
            }

            // Core size and distribution functions
            size_local_data(node_pos, node_grid);
            assign_global_to_local_data(data);
            rebase_local_data();//This could be removed with a shift in local index by halo_size in get_local_index_from_position()


            /*
            //Prints

            
            usleep(rank*1000000);

            if (rank == 0)
            {
                std::cout << "Lattice Bins [" << m_bins[0] << " " << m_bins[1] << " " << m_bins[2] << "]" << std::endl;
                std::cout << "Node Grid [" << node_grid[0] << " " << node_grid[1] << " " << node_grid[2] << "]" << std::endl;
                std::cout << "Halo [" << m_halo_size[0] << " " << m_halo_size[1] << " " << m_halo_size[2] << "]\n" << std::endl;
            }

            //std::cout << m_local_data.shape()[0] << " " << m_local_data.shape()[1] << " " << m_local_data.shape()[2] << std::endl;
            //std::cout << m_local_data.strides()[0] << " " << m_local_data.strides()[1] << " " << m_local_data.strides()[2] << std::endl;
            //std::cout << m_local_data.index_bases()[0] << " " << m_local_data.index_bases()[1] << " " << m_local_data.index_bases()[2] << std::endl;

            //Local and global bounds
            std::cout << "\nNode " << rank << " node_pos " << node_pos[0] << "/" << node_grid[0]-1 << " " << node_pos[1] << "/" << node_grid[1]-1 << " " << node_pos[2] << "/" << node_grid[2]-1 <<  std::endl;

            std::cout << "local bins  [[" << -m_halo_size[0] << ".." << m_local_block_size[0]-m_halo_size[0]-1 << "],["
                                          << -m_halo_size[1] << ".." << m_local_block_size[1]-m_halo_size[1]-1 << "],[" 
                                          << -m_halo_size[2] << ".." << m_local_block_size[2]-m_halo_size[2]-1 << "]]" << std::endl; 

            std::cout << "global bins [[" << -m_halo_size[0] + m_local_halo_offset[0] << ".." << m_local_block_size[0]-m_halo_size[0]-1 + m_local_halo_offset[0] << "],["
                                          << -m_halo_size[1] + m_local_halo_offset[1] << ".." << m_local_block_size[1]-m_halo_size[1]-1 + m_local_halo_offset[1] << "],[" 
                                          << -m_halo_size[2] + m_local_halo_offset[2] << ".." << m_local_block_size[2]-m_halo_size[2]-1 + m_local_halo_offset[2] << "]]" << std::endl; 

            //Assert test local<->global
            //double box_l[3] = {17, 3, 15.2};
            //for (double x = 0.0; x < box_l[0]; x += 0.5) {
            //    for (double y = 0.0; y < box_l[1]; y += 0.5) {
            //        for (double z = 0.0; z < box_l[2]; z += 0.5) {
            //            double p[3] = {x,y,z};

            //            T l_d;
            //            if (get_lattice_data_from_position(p,box_l,l_d))
            //                assert(data(get_global_index_from_position(p, box_l)) == l_d);
            //        }
            //    }
            //}

            //Write to file
            write_lattice_to_file();

            */
        }

        /*

        void calculate_gradient()
        {

            //shape(),index_bases() of data to boost::array
            boost::array<size_t, N> ma_shape;
            boost::array<int, N> ma_base;
            for(int dim = 0; dim < N; ++dim) {
                ma_shape[dim] = m_local_data.shape()[dim];
                ma_base[dim] = m_local_data.index_bases()[dim];
            }

            //Size gradient
            m_local_data_gradient.resize(ma_shape);
            //Rebase gradient
            m_local_data_gradient.reindex(ma_base);

            //Single loop pointer on data
            tIndexArray indexArray;

            //Loop through all elements
            for ( int i = 0; i < m_local_data.num_elements(); i++ ) {
                //The actual entry where we calculate the derivative
                indexArray = getIndexArray( m_local_data, i );
                //The gradient
                boost::array<T,N> grad;
                
                for(int dim = 0; dim < N; ++dim) {
                    //Only one bin in this dimension -> no gradient
                    if (ma_shape[dim]==1) {
                        grad[dim] = 0;
                    } else {
                        //Index arrays for neighbours
                        tIndexArray forwardIndex = indexArray;
                        tIndexArray backwardIndex = indexArray;
                        double ds = 1.0;
                        //Border: forward/backward differences
                        if (indexArray[dim] == ma_base[dim] ) {
                            forwardIndex[dim] += 1;
                        } else if (indexArray[dim] == ma_shape[dim] - 1 + ma_base[dim] ) {
                            backwardIndex[dim] -= 1;
                        //else: central differences
                        } else {
                            forwardIndex[dim] += 1;
                            backwardIndex[dim] -= 1;
                            ds = 0.5;
                        }
                        grad[dim] = ds * (m_local_data(forwardIndex)-m_local_data(backwardIndex)); //type independent minus operator needed
                    }
                }
                m_local_data_gradient(indexArray) = grad;
            }

            //write_gradient_to_file();
        }

        */
       
        /* 
        //Works for <double>
        void write_lattice_to_file()
        {
            std::ofstream file;
            file.open ("lattice.dat");
            T* p = m_local_data.data();
            tIndexArray index;
            for ( int i = 0; i < m_local_data.num_elements(); i++ ) {
                index = getIndexArray( m_local_data, i );
                //std::cout << index[0] << " " << index[1] << " " << index[2] << std::endl;
                //Indices first
                for(int dim = 0; dim < N; ++dim) {
                    file << index[dim] << " ";
                }
                //Data
                file << m_local_data(index) << std::endl;
            }
            file.close();
        }
        void write_gradient_to_file()
        {
            std::ofstream file;
            file.open ("grad.dat");
            tIndexArray index;
            for ( int i = 0; i < m_local_data_gradient.num_elements(); i++ ) {
                index = getIndexArray( m_local_data_gradient, i );
                //Indices first
                for(int dim = 0; dim < N; ++dim) {
                    file << index[dim] << " ";
                }
                //Data
                for(int dim = 0; dim < N; ++dim) {
                    file << m_local_data_gradient(index)[dim] << " ";
                }
                file << std::endl;
            }

            file.close();
        }

        */

        bool get_lattice_data_from_position(double *pos, T &data)
        {
            tIndexArray l_i;
            if (get_local_index_from_position(pos, l_i))
            {
                data = m_local_data(l_i);
                return true;
            }
            return false;
        }

        /*
        bool get_gradient_data_from_position(double *pos, boost::array<T,N> &data)
        {
            tIndexArray l_i;
            if (get_local_index_from_position(pos, l_i))
            {
                data = m_local_data_gradient(l_i);
                return true;
            }
            return false;
        }
        */



};

