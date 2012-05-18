/*
 * VTKOutputGrid.cpp
 *
 *  Created on: 23.03.2012
 *      Author: marscher
 */
#include "VTKOutputGrid.h"
using namespace ug;

bool SaveGridToVTK(Grid& grid, SubsetHandler& sh, const char* filename,
		Grid::VertexAttachmentAccessor<APosition>& aPos, int step) {
	// set max dimension of subset
	UpdateMaxDimensionOfSubset(sh, "dim");
	VTKOutputGrid output;
	return output.print_subset(&grid, &sh, filename, -1, step, 0, true);
}

template<class TElem>
void VTKOutputGrid::write_cell_types(VTKFileWriter& File, int si) {
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
//	get object type
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;

//	type
	char type;

//	get type, based on reference element type
	switch (refID) {
	case ROID_EDGE:
		type = (char) 3;
		break;
	case ROID_TRIANGLE:
		type = (char) 5;
		break;
	case ROID_QUADRILATERAL:
		type = (char) 9;
		break;
	case ROID_TETRAHEDRON:
		type = (char) 10;
		break;
	case ROID_PYRAMID:
		type = (char) 14;
		break;
	case ROID_PRISM:
		type = (char) 13;
		break;
	case ROID_HEXAHEDRON:
		type = (char) 12;
		break;
	default:
		throw(UGError("Element Type not known."));
	}

//	BStream.size = sizeof(char);

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si + 1;
	if (si < 0) {
		sistart = 0;
		siend = m_pSubsetHandler->num_subsets();
	}
	for (int si = sistart; si < siend; ++si) {
		//	get iterators
		iterBegin = m_pSubsetHandler->begin<TElem>(si);
		iterEnd = m_pSubsetHandler->end<TElem>(si);

		//	loop all elements, write type for each element to stream
		for (; iterBegin != iterEnd; ++iterBegin)
//			File.write_base64(&type);;
			File.write(&type);
	}

//	File.flush_base64_buffered();
	File.flush_base64_buffered();
}

template<class TElem>
void VTKOutputGrid::write_cell_offsets(VTKFileWriter& File, int si, int& n) {
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si + 1;
	if (si < 0) {
		sistart = 0;
		siend = m_pSubsetHandler->num_subsets();
	}
	for (int si = sistart; si < siend; ++si) {
		//	get iterators
		iterBegin = m_pSubsetHandler->begin<TElem>(si);
		iterEnd = m_pSubsetHandler->end<TElem>(si);

		//	loop all elements
		for (; iterBegin != iterEnd; ++iterBegin) {
			//	increase counter of vertices
			n += ref_elem_type::num_corners;

			//	write offset
			File.write_base64<int>(n);
		}
	}

//	flush stream
	File.flush_base64_buffered();
}

template<class TElem>
void VTKOutputGrid::write_cell_connectivity(VTKFileWriter& File, int si) {
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	get reference element id
	static const ReferenceObjectID refID = ref_elem_type::REFERENCE_OBJECT_ID;

	//	check id accessor
	UG_ASSERT(m_aaDOFIndexVRT.valid(), "ID access invalid");

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si + 1;
	if (si < 0) {
		sistart = 0;
		siend = m_pSubsetHandler->num_subsets();
	}
	for (int si = sistart; si < siend; ++si) {
		//	get iterators
		iterBegin = m_pSubsetHandler->begin<TElem>(si);
		iterEnd = m_pSubsetHandler->end<TElem>(si);

		//	loop all elements
		for (; iterBegin != iterEnd; iterBegin++) {
			//	get element
			TElem* elem = *iterBegin;

			//	write ids of the element
			if (refID != ROID_PRISM) {
				for (size_t i = 0; i < (size_t) ref_elem_type::num_corners;
						i++) {
					VertexBase* vert = elem->vertex(i);
					int id = m_aaDOFIndexVRT[vert];
					File.write_base64<int>(id);;
				}
			} else {
				int id = m_aaDOFIndexVRT[elem->vertex(0)];
				File.write_base64<int>(id);;
				id = m_aaDOFIndexVRT[elem->vertex(2)];
				File.write_base64<int>(id);;
				id = m_aaDOFIndexVRT[elem->vertex(1)];
				File.write_base64<int>(id);;
				id = m_aaDOFIndexVRT[elem->vertex(3)];
				File.write_base64<int>(id);;
				id = m_aaDOFIndexVRT[elem->vertex(5)];
				File.write_base64<int>(id);;
				id = m_aaDOFIndexVRT[elem->vertex(4)];
				File.write_base64<int>(id);;
			}
		}
	}
}

bool VTKOutputGrid::write_cells(VTKFileWriter& File, int si, int dim, int numElem,
		int numConn) {
//	write opening tag to indicate that elements will be written
	File.write( "      <Cells>\n");

	///////////////////////////
	// connectivity
	///////////////////////////
//	write opening tag to indicate that connections will be written
	File.write(
			"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	int n = sizeof(int) * numConn;
	File.write_base64<int>(n);
	File.flush_base64_buffered();

//	switch dimension
	if (numConn > 0) {
		switch (dim) {
		case 0:
			break; // no elements -> nothing to do
		case 1:
			write_cell_connectivity<Edge>(File, si);
			break;
		case 2:
			write_cell_connectivity<Triangle>(File, si);
			write_cell_connectivity<Quadrilateral>(File, si);
			break;
		case 3:
			write_cell_connectivity<Tetrahedron>(File, si);
			write_cell_connectivity<Pyramid>(File, si);
			write_cell_connectivity<Prism>(File, si);
			write_cell_connectivity<Hexahedron>(File, si);
			break;
		default:
			UG_LOG(
					"ERROR in 'VTK::write_elements': Dimension " << dim << " is not supported.\n");
			return false;
		}
	}
	File.flush_base64_buffered();

//	write closing tag
	File.write( "\n        </DataArray>\n");

	///////////////////////////
	// offsets
	///////////////////////////
//	write opening tag indicating that offsets are going to be written
	File.write(
			"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	n = sizeof(int) * numElem;
	File.write_base64<int>(n);
	File.flush_base64_buffered();

	n = 0;
//	switch dimension
	if (numElem > 0) {
		switch (dim) {
		case 0:
			break; // no elements -> nothing to do
		case 1:
			write_cell_offsets<Edge>(File, si, n);
			break;
		case 2:
			write_cell_offsets<Triangle>(File, si, n);
			write_cell_offsets<Quadrilateral>(File, si, n);
			break;
		case 3:
			write_cell_offsets<Tetrahedron>(File, si, n);
			write_cell_offsets<Pyramid>(File, si, n);
			write_cell_offsets<Prism>(File, si, n);
			write_cell_offsets<Hexahedron>(File, si, n);
			break;
		default:
			UG_LOG(
					"ERROR in 'VTK::write_elements': Dimension " << dim << " is not supported.\n");
			return false;
		}
	}
	File.flush_base64_buffered();
	File.write( "\n        </DataArray>\n");

	///////////////////////////
	// types of elements
	///////////////////////////
//	write opening tag to indicate that types will be written
	File.write(
			"        <DataArray type=\"Int8\" Name=\"types\" format=\"binary\">\n");
	File.write_base64<int>(numElem);
	File.flush_base64_buffered();

//	switch dimension
	if (numElem > 0) {
		switch (dim) {
		case 0:
			break; // no elements -> nothing to do
		case 1:
			write_cell_types<Edge>(File, si);
			break;
		case 2:
			write_cell_types<Triangle>(File, si);
			write_cell_types<Quadrilateral>(File, si);
			break;
		case 3:
			write_cell_types<Tetrahedron>(File, si);
			write_cell_types<Pyramid>(File, si);
			write_cell_types<Prism>(File, si);
			write_cell_types<Hexahedron>(File, si);
			break;
		default:
			UG_LOG(
					"ERROR in 'VTK::write_elements': Dimension " << dim << " is not supported.\n");
			return false;
		}
	}

//	write closing tag
	File.write( "\n        </DataArray>\n");
	File.write( "      </Cells>\n");

//	we're done
	return true;
}

template<typename TElem>
void VTKOutputGrid::count_sizes(int si, int& numVert, int& numElem,
		int& numConn) {
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	iterator for the elements
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	reset all marks
	m_pGrid->begin_marking();

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si + 1;
	if (si < 0) {
		sistart = 0;
		siend = m_pSubsetHandler->num_subsets();
	}
	for (int si = sistart; si < siend; ++si) {
		//	get iterators
		iterBegin = m_pSubsetHandler->begin<TElem>(si);
		iterEnd = m_pSubsetHandler->end<TElem>(si);

		//	loop over elements of this subset
		for (; iterBegin != iterEnd; ++iterBegin) {
			//	get the element
			TElem *elem = *iterBegin;

			//	count number of elements and number of connections
			++numElem;
			numConn += ref_elem_type::num_corners;

			//	loop vertices of the element
			for (size_t i = 0; i < (size_t) ref_elem_type::num_corners; ++i) {
				//	get vertex of the element
				VertexBase* v = GetVertex(elem, i);

				//	if this vertex has already been counted, skip it
				if (m_pGrid->is_marked(v))
					continue;

				// count vertex and mark it
				++numVert;
				m_pGrid->mark(v);
			}
		}
	}

//	signal end of marking
	m_pGrid->end_marking();
}
;

bool VTKOutputGrid::count_piece_sizes(int si, int dim, int& numVert,
		int& numElem, int& numConn) {
//	debug output
	UG_DLOG(LIB_DISC_OUTPUT, 2, "\n ---- Init Numbers ----\n");

//	switch dimension
	switch (dim) {
	case 0:
		if (si >= 0)
			count_sizes<VertexBase>(si, numVert, numElem, numConn);
		else {
			numVert = 0;
			for (si = 0; si < m_pSubsetHandler->num_subsets(); ++si) {
				int numVertSi = 0;
				count_sizes<VertexBase>(si, numVertSi, numElem, numConn);
				numVert += numVertSi;
			}
		}
		break;
	case 1:
		count_sizes<Edge>(si, numVert, numElem, numConn);
		break;
	case 2:
		count_sizes<Triangle>(si, numVert, numElem, numConn);
		count_sizes<Quadrilateral>(si, numVert, numElem, numConn);
		break;
	case 3:
		count_sizes<Tetrahedron>(si, numVert, numElem, numConn);
		count_sizes<Pyramid>(si, numVert, numElem, numConn);
		count_sizes<Prism>(si, numVert, numElem, numConn);
		count_sizes<Hexahedron>(si, numVert, numElem, numConn);
		break;
	default:
		UG_LOG(
				"ERROR in 'VTK::init_subset': Dimension " << dim << " is not supported.\n");
		return false;
	}

//	debug output
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Vertices: " << numVert << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Elements: " << numElem << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, "Number of Connections: " << numConn << "\n");
	UG_DLOG(LIB_DISC_OUTPUT, 2, " ---- End ----\n");

//	we're done
	return true;
}

template<typename TElem>
void VTKOutputGrid::write_points_elementwise(VTKFileWriter& File,
		Grid::VertexAttachmentAccessor<APosition>& aaPos, int si, int& n) {
//	get reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	get domain type
//	typedef typename function_type::domain_type domain_type;

//	get position attachment
//	typename domain_type::position_accessor_type& aaPos =
//			u.domain()->position_accessor();

// 	write points and remember numbering
	UG_ASSERT(m_aaDOFIndexVRT.valid(), "Missing attachment");

//	position vector
//	typename domain_type::position_type Pos;
//	typename APosition Pos;
//	typename vector3 Pos;
	vector3 Pos;

//	corner counter
	float co;

//	get iterators
//	typename TDiscreteFunction::template traits<TElem>::const_iterator iterBegin, iterEnd;
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd;

//	start marking of vertices
	m_pGrid->begin_marking();

//	loop all subsets for whole domain, only the given subset if si >= 0
	int sistart = si, siend = si + 1;
	if (si < 0) {
		sistart = 0;
		siend = m_pSubsetHandler->num_subsets(); //u.num_subsets();
	}
	for (int si = sistart; si < siend; ++si) {
		iterBegin = m_pSubsetHandler->begin<TElem>(si);
		iterEnd = m_pSubsetHandler->end<TElem>(si);

		//	loop all elements of the subset
		for (; iterBegin != iterEnd; ++iterBegin) {
			//	get the element
			TElem *elem = *iterBegin;

			//	loop vertices of the element
			for (size_t i = 0; i < (size_t) ref_elem_type::num_corners; ++i) {
				//	get vertex of element
				VertexBase* v = GetVertex(elem, i);

				//	if vertex has already be handled, skip it
				if (m_pGrid->is_marked(v))
					continue;

				//	mark the vertex as processed
				m_pGrid->mark(v);

				//	number vertex
				m_aaDOFIndexVRT[v] = n++;

				//	get position of vertex
				Pos = aaPos[v];

				//	write position to stream
				for (int i = 0; i < 3/*domain_type::dim*/; ++i) {
					co = Pos[i];
					File.write_base64<float>(co);
				}

				//	fill with missing zeros (if dim < 3)
				for (int i = 3/*domain_type::dim*/; i < 3; ++i) {
					co = 0.0;
					File.write_base64<float>(co);
				}
			}
		}
	}

	if (n > 0) {
		UG_DLOG(LIB_DISC_OUTPUT, 3,
				" ---- " << n << " Vertices (Nr. 0 - " << (n-1) << ") written to file ----\n");
	} else {
		UG_DLOG(LIB_DISC_OUTPUT, 3,
				" ---- " << n << " Vertices written to file ----\n");
	}

//	flush the stream
	File.flush_base64_buffered();

//	signal end of marking the grid
	m_pGrid->end_marking();
}

bool VTKOutputGrid::write_points(VTKFileWriter& File, int si, int dim, int numVert) {
//	write starting xml tag for points
	File.write( "      <Points>\n");
	File.write(
			"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	int n = 3 * sizeof(float) * numVert;
	File.write_base64<int>(n);
	File.flush_base64_buffered();

	Grid::VertexAttachmentAccessor<APosition> aaPos(*m_pGrid, aPosition);
//	reset counter for vertices
	n = 0;

//	switch dimension
	if (numVert > 0) {
		switch (dim) {
		case 0:
			write_points_elementwise<VertexBase>(File, aaPos, si, n);
			break;
		case 1:
			write_points_elementwise<Edge>(File, aaPos, si, n);
			break;
		case 2:
			write_points_elementwise<Triangle>(File, aaPos, si, n);
			write_points_elementwise<Quadrilateral>(File, aaPos, si, n);
			break;
		case 3:
			write_points_elementwise<Tetrahedron>(File, aaPos, si, n);
			write_points_elementwise<Pyramid>(File, aaPos, si, n);
			write_points_elementwise<Prism>(File, aaPos, si, n);
			write_points_elementwise<Hexahedron>(File, aaPos, si, n);
			break;
		default:
			UG_LOG(
					"ERROR in 'VTK::write_points': Dimension " << dim << " is not supported.\n");
			return false;
		}
	}

//	write closing tags
	File.write( "\n        </DataArray>\n");
	File.write( "      </Points>\n");

//	everything fine
	return true;
}

bool VTKOutputGrid::write_piece(VTKFileWriter& File, int si, int dim) {
//	counters
	int numVert = 0, numElem = 0, numConn = 0;

// 	Count needed sizes for vertices, elements and connections
	if (!count_piece_sizes(si, dim, numVert, numElem, numConn)) {
		UG_LOG("ERROR in 'VTK::write_subset': Can not count piece sizes.\n");
		return false;
	}

//	write the beginning of the piece, indicating the number of vertices
//	and the number of elements for this piece of the grid.
	char* str=0;
	sprintf(str, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
			numVert, numElem);
	File.write( str);

//	write vertices of this piece
	if (!write_points(File, si, dim, numVert)) {
		UG_LOG("ERROR in 'VTK::write_subset': Can not write Points.\n");
		return false;
	}

//	write elements of this piece
	if (!write_cells(File, si, dim, numElem, numConn)) {
		UG_LOG("ERROR in 'VTK::write_subset': Can not write Elements.\n");
		return false;
	}

//	write opening tag to indicate point data
	File.write( "      <PointData>\n");

//	add all components if 'selectAll' chosen
//	if (m_bSelectAll)
//	if(true)
//		for (size_t fct = 0; fct < u.num_fct(); ++fct)
//			select_nodal_scalar(u.name(fct).c_str(), u.name(fct).c_str());
//
////	loop all selected symbolic names
//	for (size_t sym = 0; sym < m_vSymbFct.size(); ++sym) {
//		//	get symb function
//		const std::string& symbNames = m_vSymbFct[sym].first;
//		const std::string& vtkName = m_vSymbFct[sym].second;
//
//		//	create function group
//		FunctionGroup fctGrp = u.fct_grp_by_name(symbNames.c_str());
//
//		//	check that all functions are contained in subset
//		bool bContained = true;
//		for (size_t i = 0; i < fctGrp.num_fct(); ++i) {
//			//	get function
//			const size_t fct = fctGrp[i];
//
//			//	check
//			if (!u.is_def_in_subset(fct, si))
//				bContained = false;
//		}
//		if (!bContained)
//			continue;

//		//	write scalar value of this function
//		if (!write_nodal_values(File, u, fctGrp, vtkName, si, dim, numVert)) {
//			UG_LOG(
//					"ERROR in 'VTK::write_subset': Can not write Scalar Values.\n");
//			return false;
//		}
//	}

//	write closing tag
	File.write( "      </PointData>\n");

//	write closing tag
	File.write( "    </Piece>\n");

//	we're done
	return true;
}

bool VTKOutputGrid::is_valid_filename(std::string& nameIn) {
// 	search for dots, they are not allowed in file name
	nameIn = nameIn.substr(0, nameIn.find_first_of('.'));

//	everything ok
	return true;
}

bool VTKOutputGrid::vtu_filename(std::string& nameOut, std::string nameIn,
		int rank, int si, int maxSi, int step) {
	if (!is_valid_filename(nameIn))
		return false;

//	copy name
	nameOut = nameIn;

//#ifdef UG_PARALLEL
//// 	process index
//	if(pcl::GetNumProcesses() > 1)
//	AppendCounterToString(nameOut, "_p", rank, pcl::GetNumProcesses() - 1);
//#endif

// 	subset index
	if (si >= 0)
		AppendCounterToString(nameOut, "_s", si, maxSi);

// 	time index
	if (step >= 0)
		AppendCounterToString(nameOut, "_t", (int) step);

// 	add file extension
	nameOut.append(".vtu");

//	we're done
	return true;
}

bool VTKOutputGrid::print_subset(Grid* grid, SubsetHandler* sh,
		const char* filename, int si, int step, number time,
		bool makeConsistent) {

	m_pGrid = grid;
	m_pSubsetHandler = sh;

//	check grid
	if (!m_pGrid) {
		UG_LOG("ERROR in 'VTK::print_subset': Cannot find underlying Grid.\n");
		return false;
	}

	//	check subsethandler
	if (!m_pSubsetHandler) {
		UG_LOG(
				"ERROR in 'VTK::print_subset': Cannot find underlying SubsetHandler.\n");
		return false;
	}

// 	attach help indices
	m_pGrid->attach_to_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.access(*m_pGrid, m_aDOFIndex);

//	get rank of process
	int rank = 0;

//	get name for *.vtu file
	std::string name;
	if (!vtu_filename(name, filename, rank, si,
			m_pSubsetHandler->num_subsets() - 1, step))
		return false;

//	open the file
	VTKFileWriter* File_p;

	try {
		File_p = new VTKFileWriter(name);
	} catch (...) {
		UG_LOG("ERROR in 'VTK::print_subset': Can not open Output File.\n");
		return false;
	}

	VTKFileWriter& File = *File_p;


//	bool if time point should be written to *.vtu file
//	in parallel we must not (!) write it to the *.vtu file, but to the *.pvtu
	bool bTimeDep = (step >= 0);
//#ifdef UG_PARALLEL
//	if(pcl::GetNumProcesses() > 1) bTimeDep = false;
//#endif
//	header
	File.write( "<?xml version=\"1.0\"?>\n");
	char* str=0;
	sprintf(str,
	 "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
			"byte_order=\"%s\">\n",
#ifdef __SWAPBYTES__
			"LittleEndian"
#else
			"BigEndian"
#endif
);
	File.write(str);

//	writing time point
	if (bTimeDep) {
		char* str=0;
		sprintf(str,"  <Time timestep=\"%g\"/>\n",time );
		File.write(str);
	}

//	opening the grid
	File.write( "  <UnstructuredGrid>\n");

// 	get dimension of grid-piece
	int dim = -1;
	if (si >= 0)
		dim = DimensionOfSubset(*m_pSubsetHandler, si);
	else
		dim = DimensionOfSubsets(*m_pSubsetHandler);

//	write piece of grid
	if (dim >= 0) {
//		if (!write_piece(File, u, si, dim)) {
		if (!write_piece(File, si, dim)) {
			UG_LOG("ERROR in 'VTK::print_subset': Can not write Subset.\n");
			return false;
		}
	} else {
		//	if dim < 0, some is wrong with grid, except no element is in the grid
		if (((si < 0) && m_pGrid->num<VertexBase>() != 0)
				|| ((si >= 0) && m_pSubsetHandler->num<VertexBase>(si) != 0)) {
			UG_LOG("ERROR in 'VTK::print_subset': Dimension of grid/subset not"
			" detected correctly although grid objects present.\n");
			return false;
		}

		//	write that no elements are in the grid
		char* str2 =0;
		sprintf(str2,
				"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 0,
				0);
		File.write(str2);
		if (!write_points(File, si, dim, 0)
				|| !write_cells(File, si, dim, 0, 0)) {
			UG_LOG(
					"ERROR in 'VTK::print_subset': Can not write empty Points/Cells.\n");
			return false;
		}
		File.write( "    </Piece>\n");

	}
//	write closing xml tags
	File.write( "  </UnstructuredGrid>\n");
	File.write( "</VTKFile>\n");

// 	detach help indices
	m_pGrid->detach_from_vertices(m_aDOFIndex);
	m_aaDOFIndexVRT.invalidate();

//#ifdef UG_PARALLEL
////	write grouping *.pvtu file in parallel case
//	if(!write_pvtu(u, filename, si, step, time))
//	{
//		UG_LOG("ERROR in 'VTK::print_subset': Can not write pvtu - file.\n");
//		fclose(File);
//		return false;
//	}
//#endif

//	remember time step
//	static std::map<std::string, std::vector<number> > m_mTimestep;
//	if (step >= 0) {
//		//	get vector of time points for the name
//		std::vector<number> &vTimestep = m_mTimestep[std::string(filename)];
//
//		//	resize the vector
//		vTimestep.resize(step + 1);
//
//		//	write time point
//		vTimestep[step] = time;
//	}

//	we're done
	delete File_p;
	return true;
}
