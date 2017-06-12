inline uint minComponentId(void) const {uint i; this->minCoeff(&i); return i;}
inline uint maxComponentId(void) const {uint i; this->maxCoeff(&i); return i;}


// 6D scalar accessors
inline Scalar at(uint idx_alpha1, uint idx_theta1, uint idx_phi1,
				 uint idx_alpha2, uint idx_theta2, uint idx_phi2) const
{
	uint idx1 = idx_alpha1*N_THETA*N_PHI + idx_theta1*N_PHI + idx_phi1;
	uint idx2 = idx_alpha2*N_THETA*N_PHI + idx_theta2*N_PHI + idx_phi2;
	
	return this->operator()(idx1, idx2);
}

inline Scalar& at(uint idx_alpha1, uint idx_theta1, uint idx_phi1,
				  uint idx_alpha2, uint idx_theta2, uint idx_phi2)
{
	uint idx1 = idx_alpha1*N_THETA*N_PHI + idx_theta1*N_PHI + idx_phi1;
	uint idx2 = idx_alpha2*N_THETA*N_PHI + idx_theta2*N_PHI + idx_phi2;
	
	return this->operator()(idx1, idx2);
}

// 3D scalar accessors
inline Scalar at(uint idx_alpha, uint idx_theta, uint idx_phi) const
{
	uint idx = idx_alpha*N_THETA*N_PHI + idx_theta*N_PHI + idx_phi;
	
	return this->operator()(idx);
}

inline Scalar& at(uint idx_alpha, uint idx_theta, uint idx_phi)
{
	uint idx = idx_alpha*N_THETA*N_PHI + idx_theta*N_PHI + idx_phi;
	
	return this->operator()(idx);
}

// 3D row accessors
inline ConstRowXpr row_at(uint idx_alpha, uint idx_theta, uint idx_phi) const
{
	uint idx = idx_alpha*N_THETA*N_PHI + idx_theta*N_PHI + idx_phi;
	
	return this->row(idx);
}

inline RowXpr row_at(uint idx_alpha, uint idx_theta, uint idx_phi)
{
	uint idx = idx_alpha*N_THETA*N_PHI + idx_theta*N_PHI + idx_phi;
	
	return this->row(idx);
}

// 3D col accessors
inline ConstColXpr col_at(uint idx_alpha, uint idx_theta, uint idx_phi) const
{
	uint idx = idx_alpha*N_THETA*N_PHI + idx_theta*N_PHI + idx_phi;
	
	return this->col(idx);
}

inline ColXpr col_at(uint idx_alpha, uint idx_theta, uint idx_phi)
{
	uint idx = idx_alpha*N_THETA*N_PHI + idx_theta*N_PHI + idx_phi;
	
	return this->col(idx);
}

// Symmetrise matrix
inline DenseBase<Derived>& symmetrise()
{
	for ( uint idx_row = 0; idx_row < E_DIM; ++idx_row )
	{
		for ( uint idx_col = 0; idx_col < idx_row; ++idx_col )
		{
			this->operator()(idx_row, idx_col) = (this->operator()(idx_row, idx_col) + this->operator()(idx_col, idx_row)) /2.;
			this->operator()(idx_col, idx_row) = this->operator()(idx_row, idx_col);
		}
	}
	
	return *this;
}
