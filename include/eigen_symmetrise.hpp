inline uint minComponentId(void) const  {int i; this->minCoeff(&i); return i;}
inline uint maxComponentId(void) const  {int i; this->maxCoeff(&i); return i;}


inline uint get_sym_size(uint idx)       {return idx*(idx+1)/2;};
inline uint get_sym_size(uint idx) const {return idx*(idx+1)/2;};


inline Scalar sym(uint idx_row, uint idx_col, uint dim) const
{
	eigen_assert( (idx_row < dim) && (idx_col < dim) );
	uint idx = (idx_row <= idx_col) ? get_sym_size(idx_col) + idx_row : get_sym_size(idx_row) + idx_col;
	
	return this->operator()(idx);
}

inline Scalar& sym(uint idx_row, uint idx_col, uint dim)
{
	eigen_assert( (idx_row < dim) && (idx_col < dim) );
	uint idx = (idx_row <= idx_col) ? get_sym_size(idx_col) + idx_row : get_sym_size(idx_row) + idx_col;
	
	return this->operator()(idx);
}

inline ArrayBase<Derived>& sym_normalise(uint dim)
{
	for ( uint idx_row = 0; idx_row < dim; ++idx_row )
	{
		for ( uint idx_col = 0; idx_col < idx_row; ++idx_col )
		{
			this->sym(idx_row, idx_col, dim) /= 2.;
		}
	}
	
	return *this;
}
