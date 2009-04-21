#include "tasks.hpp"

#include <deque>

using namespace std;

Tasks::Tasks(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode) {

	if (mode=="none") {
		generate_independent_tasks(frames,qqq,nn,mode);
	}
	else {
		generate_tasks_by_framecoupling(frames,qqq,nn,mode);
	}
};

void Tasks::print() {
		for(size_t i = 0; i < size(); ++i)
		{
			Task& t = this->at(i);
			clog << "INFO>> " << "(id=" << t.id << ") " << "q=" << t.q << "\t" << "{";
			for(size_t j = 0; j < t.rftable.size(); ++j)
			{
				pair<int,int> rf = t.rftable[j];
				clog << rf.first << ",";
			}
			clog << "}" << "=>" << "{";
			for(size_t j = 0; j < t.rftable.size(); ++j)
			{
				pair<int,int> rf = t.rftable[j];
				clog << rf.second << ",";
			}	
			clog << "}" ;
			clog << endl;
		}
}

void Tasks::generate_independent_tasks(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode)  {
	// nn = number of nodes in total

	size_t nf = frames.size();
	size_t nq = qqq.size();

	size_t current_task_id=0;

	deque<int> idlenodes;
	// make all nodes idle;
	for(size_t i = 0; i < nn; ++i)
	{
		idlenodes.push_back(i);
	}
	
	size_t q=0;
	while (q<nq) {
		
		for(size_t f = 0; f < nf; ++f)
		{
			// try to increase total idlenodes to at least number of frames
			while (idlenodes.size()<nf) {
				// put next set of idle nodes onto deque
				for(size_t i = 0; i < nn; ++i)
				{
					idlenodes.push_back(i);
				}
			}
			


			Task tm;
			tm.q = qqq[q];
			int thisnode = idlenodes.front();
			idlenodes.pop_front();
			if ( find(tm.ranks.begin(),tm.ranks.end(),thisnode) == tm.ranks.end() ) tm.ranks.push_back(thisnode);
			tm.rftable.push_back(make_pair(thisnode,f));
			tm.id = current_task_id++;
			tm.mode = mode;
			push_back(tm);
			
		}	

		q++;

// not necessary for independent tasks:
		// align nodes to frame boundary
//		if ( idlenodes.size() < nf ) {
//			idlenodes.clear();
//		}
	}	
}


void Tasks::generate_tasks_by_framecoupling(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode) {
	// nn = number of nodes in total

	size_t nf = frames.size();
	size_t nq = qqq.size();

	size_t current_task_id=0;

	deque<int> idlenodes;
	// make all nodes idle;
	for(size_t i = 0; i < nn; ++i)
	{
		idlenodes.push_back(i);
	}
	
	size_t q=0;
	while (q<nq) {
		
		// try to increase total idlenodes to at least number of frames
		while (idlenodes.size()<nf) {
			// put next set of idle nodes onto deque
			for(size_t i = 0; i < nn; ++i)
			{
				idlenodes.push_back(i);
			}
		}
			
		Task tm;
		tm.q = qqq[q];
		for(size_t f = 0; f < nf; ++f)
		{	
			int thisnode = idlenodes.front();
			idlenodes.pop_front();
			if ( find(tm.ranks.begin(),tm.ranks.end(),thisnode) == tm.ranks.end() ) tm.ranks.push_back(thisnode);
			tm.rftable.push_back(make_pair(thisnode,f));
		}	
		tm.id = current_task_id++;
		tm.mode = mode;
		push_back(tm);
		q++;

		// align nodes to frame boundary
		if ( idlenodes.size() < nf ) {
			idlenodes.clear();
		}
	}
}

Tasks Tasks::scope(int rank) {
	Tasks result;
	for(size_t i = 0; i < size(); ++i)
	{
		Task& t = this->at(i);
		for(size_t j = 0; j < t.rftable.size(); ++j)
		{
			pair<int,int> rf = t.rftable[j];
			if (rf.first == rank) {
				result.push_back(t); 
				break;
			};
		}
	}
	return result;
}