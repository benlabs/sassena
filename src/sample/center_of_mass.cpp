/** \file
This file contains a class which defines represents a center of mass calculation based on coordinate input, and related structures used for fitting.

\author Benjamin Lindner <ben@benlabs.net>
\version 1.3.0
\copyright GNU General Public License
*/
 
// direct header
#include "sample/center_of_mass.hpp"

// other headers
#include "../vendor/boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "../vendor/boost/numeric/bindings/lapack/gesvd.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "sample/atoms.hpp"
#include "sample/atomselection.hpp"
#include "math/coor3d.hpp"
#include "sample/coordinate_set.hpp"
#include "control.hpp"
#include "log.hpp"

using namespace std;

// code snippets for determinant are from http://www.anderswallin.net/

inline int determinant_sign(const boost::numeric::ublas::permutation_matrix<size_t>& pm)
{
    int pm_sign=1;
    std::size_t size = pm.size();
    for (std::size_t i = 0; i < size; ++i)
        if (i != pm(i))
            pm_sign *= -1.0; // swap_rows would swap a pair of rows here, so we change sign
    return pm_sign;
}
 
inline double determinant(boost::numeric::ublas::matrix<double>& m ) {
    boost::numeric::ublas::permutation_matrix<size_t> pm(m.size1());
    double det = 1.0;
    if( boost::numeric::ublas::lu_factorize(m,pm) ) {
        det = 0.0;
    } else {
        for(size_t i = 0; i < m.size1(); i++) 
            det *= m(i,i); // multiply by elements on diagonal
        det = det * determinant_sign( pm );
    }
    return det;
}

Fit::Fit(
    Atoms& atoms, // contains IDs in sequence
    CoordinateSet& cs, // contains original coordinates of target
    IAtomselection* pcs_selection, // selection corresponding to target
    IAtomselection* pcs_selection_manip, // (sub-)selection of target which specifies which atoms to move
    CoordinateSet& cs_ref,  // coordinate set of the reference structure
    IAtomselection* pref_selection ) // selection correspoding to reference
{
    IAtomselection& cs_selection_manip = *pcs_selection_manip;
    IAtomselection& ref_selection = *pref_selection;
    
    // since rotational fitting has some algorithmic complexity, we
    // create a reduced copy for the target coordinate set to work with
    // this makes things easier, i.e. we can neglect some mapping

    CoordinateSet cs_redcopy(cs,pcs_selection,pref_selection);
    if (cs_redcopy.size()!=cs_ref.size()) {
        Err::Inst()->write("Fitting requires the reference and the target to contain the same number of atoms");
        throw;
    }
    
    // pre-determine center of mass vectors for each coordinate set
    CartesianCoor3D ref = CenterOfMass(atoms,cs_ref,pref_selection);
    CartesianCoor3D pos = CenterOfMass(atoms,cs_redcopy,pref_selection);
    
    // compute the correlation kernel
    using namespace boost::numeric::ublas;
    
    matrix<double> K = zero_matrix<double>(3, 3);
    boost::numeric::ublas::vector<double> v(3);
    boost::numeric::ublas::vector<double> u(3);
    for(size_t i = 0; i < ref_selection.size(); ++i)
    {
        size_t index = ref_selection[i];
        size_t id = atoms[index];
        double mi = Database::Inst()->masses.get(id);
        v[0]=cs_redcopy.c1[i] - pos.x;
        v[1]=cs_redcopy.c2[i] - pos.y;
        v[2]=cs_redcopy.c3[i] - pos.z;
        u[0]=cs_ref.c1[i] - ref.x;
        u[1]=cs_ref.c2[i] - ref.y;
        u[2]=cs_ref.c3[i] - ref.z;
//        cout << "ref, " << cs_ref.c1[i] << ", "  << cs_ref.c2[i] << ", "  << cs_ref.c3[i] << endl;
        
        K +=mi*outer_prod(u, v);
    }
    boost::numeric::ublas::matrix<double> U(3,3);
    boost::numeric::ublas::matrix<double> Vt(3,3);    
    boost::numeric::ublas::vector<double> S(3);

    // do singular value decomposition for K
//    boost::numeric::bindings::lapack::gesvd(K,S,U,Vt);
    // U and Vt are swapped!
    boost::numeric::bindings::lapack::gesvd(K,S,Vt,U);
    
    // now the crucial information is contained in the determinant of the dot product between U and Vt
    boost::numeric::ublas::matrix<double> UVt = prod(U,Vt);
    double detUVt = determinant(UVt);
    boost::numeric::ublas::matrix<double> P = zero_matrix<double>(3,3);
    P(0,0)=1;P(1,1)=1;P(2,2)=detUVt;
    boost::numeric::ublas::matrix<double> PVt = prod(P,Vt);
    boost::numeric::ublas::matrix<double> UPVt = prod( U , PVt );
    
    // UPVt is the rotation matrix, rotate & translate
    boost::numeric::ublas::vector<double> before(3);
    boost::numeric::ublas::vector<double> after(3);
    for(size_t i = 0; i < ref_selection.size(); ++i)
    {
        before[0] = cs_redcopy.c1[i] - pos.x;
        before[1] = cs_redcopy.c2[i] - pos.y;
        before[2] = cs_redcopy.c3[i] - pos.z;
        
//        cout << "in, " << cs_redcopy.c1[i] << ", "  << cs_redcopy.c2[i] << ", "  << cs_redcopy.c3[i] << endl;
        
        after = prod(UPVt,before);
        cs_redcopy.c1[i]=after[0]+ref.x;
        cs_redcopy.c2[i]=after[1]+ref.y;
        cs_redcopy.c3[i]=after[2]+ref.z;
    }
    
    // now the coordinates are in the local copy, we need  to map them back into the original coordinate set.
    size_t csel_total = cs_selection_manip.size();
    size_t ssel_total = ref_selection.size();
    size_t csel_iter = 0; // iterates target coordinate set
    size_t ssel_iter = 0; // iterates local copy

    while( ( csel_iter < csel_total) && (ssel_iter < ssel_total) ) {
        size_t csel_index = cs_selection_manip[csel_iter];
        size_t ssel_index = ref_selection[ssel_iter];

        if (csel_index==ssel_index) {
            cs.c1[csel_index]=cs_redcopy.c1[ssel_index];
            cs.c2[csel_index]=cs_redcopy.c2[ssel_index];
            cs.c3[csel_index]=cs_redcopy.c3[ssel_index];
//            cout << "out, " << cs.c1[csel_index] << ", "  << cs.c2[csel_index] << ", "  << cs.c3[csel_index] << endl;
            ssel_iter++;
            csel_iter++;           
        } else if (csel_index>ssel_index) {
            ssel_iter++;
        } else if (csel_index<ssel_index) {
            csel_iter++;
        }
    }
    //done
}

// Determines the COFM for a sub-selection of a coordinate set (which may itself be a subselection of the entire system)
CenterOfMass::CenterOfMass(Atoms& atoms,CoordinateSet& cs, IAtomselection* pcs_selection,  IAtomselection* pcofm_selection) {

    IAtomselection& cs_selection = *pcs_selection;
    IAtomselection& cofm_selection = *pcofm_selection;
    
    if (cs.get_representation()!=CARTESIAN) {
        Err::Inst()->write("Center of Mass only implementated for cartesian Coordinate Sets");
        throw;
    }

    size_t csel_total = cs_selection.size();
    size_t ssel_total = cofm_selection.size();
    
        if (ssel_total<1) {
            Warn::Inst()->write("Warning! Computing Center of Mass for an empty atomselection");
            Warn::Inst()->write("Setting Center of mass to (0,0,0)");
            m_center = CartesianCoor3D(0,0,0);
        } else if (csel_total<1) {
            Warn::Inst()->write("Warning! Computing Center of Mass for an empty coordinate set");
            Warn::Inst()->write("Setting Center of mass to (0,0,0)");
            m_center = CartesianCoor3D(0,0,0);          
    } else {

        coor2_t xt,yt,zt;
        double m,mi;
        xt = yt = zt = 0.0;
        m = 0.0;
        
        size_t csel_iter = 0;
        size_t ssel_iter = 0;
        
        while( ( csel_iter < csel_total) && (ssel_iter < ssel_total) ) {
            size_t csel_index = cs_selection[csel_iter];
            size_t ssel_index = cofm_selection[ssel_iter];
        
            if (csel_index==ssel_index) {
                size_t csel_id = atoms[csel_index];
                mi = Database::Inst()->masses.get(csel_id);
                        m += mi;
        
                        xt += cs.c1[csel_iter]*mi;
                        yt += cs.c2[csel_iter]*mi;
                        zt += cs.c3[csel_iter]*mi;
        
                ssel_iter++;
                csel_iter++;           
            } else if (csel_index>ssel_index) {
                ssel_iter++;
            } else if (csel_index<ssel_index) {
                csel_iter++;
            }
        }
        
        m_center = CartesianCoor3D(xt/m,yt/m,zt/m);
        
    }

}

CenterOfMass::CenterOfMass(Atoms& atoms,Frame& frame,IAtomselection* pselection) {
	
    IAtomselection& selection = *pselection;
	
	if (selection.size()==0) {
		Err::Inst()->write("Warning! Computing Center of Mass for an empty atomselection");
		Err::Inst()->write("Setting Center of mass to (0,0,0)");		
		m_center = CartesianCoor3D(0,0,0);
        return;
	}

    coor2_t xt,yt,zt;
    double m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	for (size_t i=0;i<selection.size();i++) {
        size_t sel_id = atoms[selection[i]];
        mi = Database::Inst()->masses.get(sel_id);
		m += mi;
		xt += frame.x[selection[i]]*mi;
		yt += frame.y[selection[i]]*mi;
		zt += frame.z[selection[i]]*mi;
	}
	
	m_center =  CartesianCoor3D(xt/m,yt/m,zt/m);
}

// Determines the COFM for the entire selection of a coordinate set
CenterOfMass::CenterOfMass(Atoms& atoms,CoordinateSet& cset,IAtomselection* pselection) {
	
    IAtomselection& selection = *pselection;
	
	if (selection.size()==0) {
		Err::Inst()->write("Warning! Computing Center of Mass for an empty atomselection");
		Err::Inst()->write("Setting Center of mass to (0,0,0)");		
		m_center = CartesianCoor3D(0,0,0);
        return;
	}

    coor2_t xt,yt,zt;
    double m,mi;
	xt = yt = zt = 0.0;
	m = 0.0;

	for (size_t i=0;i<selection.size();i++) {
        size_t sel_id = atoms[selection[i]];
        mi = Database::Inst()->masses.get(sel_id);
		m += mi;
		xt += cset.c1[i]*mi;
		yt += cset.c2[i]*mi;
		zt += cset.c3[i]*mi;
	}
	
	m_center =  CartesianCoor3D(xt/m,yt/m,zt/m);
}

CenterOfMass::operator CartesianCoor3D() {
	return m_center;
}

// end of file
