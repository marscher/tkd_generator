/*
 * VTKOutputGrid.h
 *
 *  Created on: 21.03.2012
 *      Author: marscher
 */
#ifndef __vtk_output_grid__
#define __vtk_output_grid__
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_grid/lib_grid.h"
#include "lib_disc/io/vtkoutput.h"
#include "ugbase.h"

using namespace ug;

bool SaveGridToVTK(Grid& grid, SubsetHandler& sh, const char* filename,
		Grid::VertexAttachmentAccessor<APosition>& aPos, int step);

class VTKOutputGrid {

public:

	bool print_subset(Grid* grid, SubsetHandler* sh, const char* filename,
			int si, int step, number time, bool makeConsistent);

	///	default constructor
	VTKOutputGrid() :
			m_bSelectAll(true), m_pGrid(NULL), m_pSubsetHandler(NULL) {
	}

protected:
	template <typename TElem> void
	write_cell_types(FILE* File, int si);

	template <typename TElem>void
	write_cell_offsets(FILE* File, int si, int& n);

	template <typename TElem>void
	write_cell_connectivity(FILE* File, int si);

	bool write_cells(FILE* File, int si, int dim, int numElem, int numConn);

	template <typename TElem> void count_sizes(int si, int& numVert, int& numElem, int& numConn);

	bool count_piece_sizes(int si, int dim, int& numVert, int& numElem,
			int& numConn);

	template <typename TElem> void write_points_elementwise(FILE* File,
			Grid::VertexAttachmentAccessor<APosition>& aaPos, int si, int& n);

	bool write_points(FILE* File, int si, int dim, int numVert);

	bool write_piece(FILE* File, int si, int dim);

	bool is_valid_filename(std::string& nameIn);

	bool vtu_filename(std::string& nameOut, std::string nameIn, int rank,
			int si, int maxSi, int step);

///	scheduled components to be printed
	bool m_bSelectAll;

///	Grid, working on
	Grid* m_pGrid;
	SubsetHandler* m_pSubsetHandler;

///	attachment for dofs
	typedef ug::Attachment<int> ADOFIndex;
	ADOFIndex m_aDOFIndex;
	Grid::VertexAttachmentAccessor<ADOFIndex> m_aaDOFIndexVRT;
};
#endif
