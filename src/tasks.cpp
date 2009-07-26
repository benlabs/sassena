#include "tasks.hpp"

#include <deque>

using namespace std;

int Task::timeunits() {
	map<int,int> nodecounttable;
	for(size_t i = 0; i < rftable.size(); ++i)
	{
		if (nodecounttable.find(rftable[i].first)!=nodecounttable.end()) {
			nodecounttable[rftable[i].first] = nodecounttable[rftable[i].first] + 1;
		}
		else {
			nodecounttable[rftable[i].first] = 1;
		}
	}
	int max=1;
	for(map<int,int>::iterator ncti = nodecounttable.begin(); ncti!=nodecounttable.end();ncti++)
	{
		max = (ncti->second>max) ? ncti->second : max;
	}
	return max;
}

std::vector<int> Task::frames(int rank) {
	vector<int> result;
	for(size_t i = 0; i < rftable.size(); ++i)
	{
		if (rftable[i].first==rank) result.push_back(rftable[i].second);
	}
	return result;
}

size_t Task::frames_max() {

	map<int,size_t> count;
	for(size_t i = 0; i < rftable.size(); ++i)
	{
		count[rftable[i].first]+=1;
	}
	size_t max=0;
	for(map<int,size_t>::iterator ci = count.begin(); ci!=count.end(); ++ci)
	{
		if (max<ci->second) max=ci->second;
	}
	return max;
}


Tasks::Tasks(std::vector<int>& frames,std::vector<CartesianCoor3D>& qqq,size_t nn,std::string mode) {

	if (mode=="time") {
		generate_tasks_by_framecoupling(frames,qqq,nn,mode);
	}
	else if (mode=="none") {
		generate_independent_tasks(frames,qqq,nn,mode);
	}
	else {
		cerr << "ERROR>> " << "Correlation type not understood" << endl;
		throw;
	}
};

void Tasks::print() {
//		for(size_t i = 0; i < size(); ++i)
//		{
//			Task& t = this->at(i);
//			clog << "INFO>> " << "(id=" << t.id << "," << "color=" << t.color << ") " << "q=" << t.q << "\t" << "{";
//			for(size_t j = 0; j < t.rftable.size(); ++j)
//			{
//				pair<int,int> rf = t.rftable[j];
//				clog << rf.first << ",";
//			}
//			clog << "}" << "=>" << "{";
//			for(size_t j = 0; j < t.rftable.size(); ++j)
//			{
//				pair<int,int> rf = t.rftable[j];
//				clog << rf.second << ",";
//			}	
//			clog << "}" ;
//			clog << endl;
//		}
		clog << "INFO>> " << "Total number of tasks: " << size() << endl;
		clog << "INFO>> " << "Number of simultaneous jobs allowed: " << colors << endl;
		
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
	size_t color =0;
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
			tm.color = 0;
			push_back(tm);
			
		}	

		q++; 

	}
	colors = nn; // indicates the number of simultaneous indepentend tasks	
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
	size_t color=0;
	size_t max_color = 0;
	while (true) {
		
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
		tm.color = color;
		push_back(tm);
		q++; if (q>=nq) break;
		color++;

		// align nodes to frame boundary
		if ( idlenodes.size() < nf ) {
			idlenodes.clear();
			if (max_color<color) max_color=color;
			color=0;			
		}
	}
	colors = max_color;
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