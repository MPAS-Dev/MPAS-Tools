
class edge {/*{{{*/
	public:
		int cell1, cell2;
		int vertex1, vertex2;
		int idx;

	edge() 
		: cell1(-1), cell2(-1), vertex1(-1), vertex2(-1), idx(-1) { }

	struct edge_hasher {/*{{{*/
		size_t operator()(const edge &edge_) const {
			uint32_t hash; 
			size_t i, key[2] = { (size_t)edge_.vertex1, (size_t)edge_.vertex2 };
			for(hash = i = 0; i < sizeof(key); ++i) {
				hash += ((uint8_t *)key)[i];
				hash += (hash << 10);
				hash ^= (hash >> 6);
			}
			hash += (hash << 3);
			hash ^= (hash >> 11);
			hash += (hash << 15);
			return hash;
		}
	};/*}}}*/

	bool operator==(const edge &e) const {/*{{{*/
		return (vertex1 == e.vertex1) && (vertex2 == e.vertex2);
	}/*}}}*/
};/*}}}*/


