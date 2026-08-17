#include "TabataAssemblyGPU.h"
#include <cuda_runtime.h>
using HYPRE_BigInt = int; using HYPRE_Real = double;

__global__ void TabataAssemblyKernel(int numTet, int dim, int local_ndofs, int my_rank,
    const int* elem_dofs, const double* dshape, const double* volume,
    const HYPRE_BigInt* global_dofs, const int* isBC, const double* rho_prev,
    const double* velocity, const HYPRE_BigInt* upwind_gdofs,
    const double* upwind_dshape, const int* row_ptr, const HYPRE_BigInt* col_idx,
    double* vals) {
    int tid=blockIdx.x*blockDim.x+threadIdx.x;
    if(tid>=numTet*4)return;
    int e=tid/4,i=tid%4;
    int ld=elem_dofs[e*4+i];
    HYPRE_BigInt gd=global_dofs[ld];
    int r0=row_ptr[ld],r1=row_ptr[ld+1];
    if(isBC[ld]){for(int k=r0;k<r1;k++)if(col_idx[k]==gd){vals[k]=1.0;return;}return;}
    double mass=rho_prev[ld]*volume[e]/4.0;
    for(int k=r0;k<r1;k++)if(col_idx[k]==gd){atomicAdd(&vals[k],mass);break;}
    int ub=ld*4;
    if(upwind_gdofs[ub]<0)return;
    for(int j=0;j<4;j++){HYPRE_BigInt gj=upwind_gdofs[ub+j];double conv=0.0;
        for(int d=0;d<dim;d++)conv+=velocity[ld*dim+d]*upwind_dshape[(ub+j)*dim+d];
        double val=conv*volume[e]/4.0;
        for(int k=r0;k<r1;k++)if(col_idx[k]==gj){atomicAdd(&vals[k],val);break;}}
}

// ── 静态缓存：mesh 拓扑相关的 GPU buffer（跨迭代复用） ──
static int    *d_ed_static  = nullptr; static int    ed_static_size  = 0;
static double *d_ds_static  = nullptr; static int    ds_static_size  = 0;
static double *d_vol_static = nullptr; static int    vol_static_size = 0;
static GIdx   *d_gd_static  = nullptr; static int    gd_static_size  = 0;
static int    *d_bc_static  = nullptr; static int    bc_static_size  = 0;

// ── 静态缓存：每次迭代值变但大小不变的 buffer ──
static double *d_rp_buf   = nullptr; static int    rp_buf_size   = 0;
static double *d_vl_buf   = nullptr; static int    vl_buf_size   = 0;
static GIdx   *d_ug_buf   = nullptr; static int    ug_buf_size   = 0;
static double *d_ud_buf   = nullptr; static int    ud_buf_size   = 0;
static int    *d_rptr_buf = nullptr; static int    rptr_buf_size = 0;
static GIdx   *d_ci_buf   = nullptr; static int    ci_buf_size   = 0;

static void ensure_capacity(void **dp, int *size, int need_bytes) {
    if (*dp == nullptr || *size < need_bytes) {
        if (*dp) cudaFree(*dp);
        cudaMalloc(dp, need_bytes);
        *size = need_bytes;
    }
}

void TabataGPUComputeValues(int numTet,int dim,int local_ndofs,int my_rank,
    const std::vector<int>&h_ed, const std::vector<double>&h_ds, const std::vector<double>&h_vol,
    const std::vector<GIdx>&h_ug, const std::vector<double>&h_ud,
    const std::vector<GIdx>&h_gd, const std::vector<int>&h_bc,
    const std::vector<double>&h_rp, const std::vector<double>&h_vl,
    const std::vector<int>&h_rptr, const std::vector<GIdx>&h_ci,
    double**d_vals_out,int*nnz_out){
    int nnz=(int)h_ci.size();*nnz_out=nnz;
    double*d_vals;cudaMalloc(&d_vals,nnz*sizeof(double));cudaMemset(d_vals,0,nnz*sizeof(double));

    // ── 不变数据：首次或大小变化时分配 + H2D ──
    int ed_bytes = numTet*4*sizeof(int);
    if (d_ed_static == nullptr || ed_static_size < ed_bytes) {
        if (d_ed_static) cudaFree(d_ed_static);
        cudaMalloc(&d_ed_static, ed_bytes);
        cudaMemcpy(d_ed_static, h_ed.data(), ed_bytes, cudaMemcpyHostToDevice);
        ed_static_size = ed_bytes;
    }

    int ds_bytes = numTet*4*dim*sizeof(double);
    if (d_ds_static == nullptr || ds_static_size < ds_bytes) {
        if (d_ds_static) cudaFree(d_ds_static);
        cudaMalloc(&d_ds_static, ds_bytes);
        cudaMemcpy(d_ds_static, h_ds.data(), ds_bytes, cudaMemcpyHostToDevice);
        ds_static_size = ds_bytes;
    }

    int vol_bytes = numTet*sizeof(double);
    if (d_vol_static == nullptr || vol_static_size < vol_bytes) {
        if (d_vol_static) cudaFree(d_vol_static);
        cudaMalloc(&d_vol_static, vol_bytes);
        cudaMemcpy(d_vol_static, h_vol.data(), vol_bytes, cudaMemcpyHostToDevice);
        vol_static_size = vol_bytes;
    }

    int gd_bytes = local_ndofs*sizeof(GIdx);
    if (d_gd_static == nullptr || gd_static_size < gd_bytes) {
        if (d_gd_static) cudaFree(d_gd_static);
        cudaMalloc(&d_gd_static, gd_bytes);
        cudaMemcpy(d_gd_static, h_gd.data(), gd_bytes, cudaMemcpyHostToDevice);
        gd_static_size = gd_bytes;
    }

    int bc_bytes = local_ndofs*sizeof(int);
    if (d_bc_static == nullptr || bc_static_size < bc_bytes) {
        if (d_bc_static) cudaFree(d_bc_static);
        cudaMalloc(&d_bc_static, bc_bytes);
        cudaMemcpy(d_bc_static, h_bc.data(), bc_bytes, cudaMemcpyHostToDevice);
        bc_static_size = bc_bytes;
    }

    // ── 可变数据：复用 buffer，只做 H2D ──
    ensure_capacity((void**)&d_rp_buf,   &rp_buf_size,   local_ndofs*sizeof(double));
    ensure_capacity((void**)&d_vl_buf,   &vl_buf_size,   local_ndofs*dim*sizeof(double));
    ensure_capacity((void**)&d_ug_buf,   &ug_buf_size,   local_ndofs*4*sizeof(GIdx));
    ensure_capacity((void**)&d_ud_buf,   &ud_buf_size,   local_ndofs*4*dim*sizeof(double));
    ensure_capacity((void**)&d_rptr_buf, &rptr_buf_size,  (local_ndofs+1)*sizeof(int));
    ensure_capacity((void**)&d_ci_buf,   &ci_buf_size,    nnz*sizeof(GIdx));

    cudaMemcpy(d_rp_buf,   h_rp.data(),   local_ndofs*sizeof(double),         cudaMemcpyHostToDevice);
    cudaMemcpy(d_vl_buf,   h_vl.data(),   local_ndofs*dim*sizeof(double),     cudaMemcpyHostToDevice);
    cudaMemcpy(d_ug_buf,   h_ug.data(),   local_ndofs*4*sizeof(GIdx),         cudaMemcpyHostToDevice);
    cudaMemcpy(d_ud_buf,   h_ud.data(),   local_ndofs*4*dim*sizeof(double),   cudaMemcpyHostToDevice);
    cudaMemcpy(d_rptr_buf, h_rptr.data(), (local_ndofs+1)*sizeof(int),        cudaMemcpyHostToDevice);
    cudaMemcpy(d_ci_buf,   h_ci.data(),   nnz*sizeof(GIdx),                   cudaMemcpyHostToDevice);

    int threads=256,blocks=(numTet*4+255)/256;
    TabataAssemblyKernel<<<blocks,threads>>>(numTet,dim,local_ndofs,my_rank,
        d_ed_static,d_ds_static,d_vol_static,d_gd_static,d_bc_static,
        d_rp_buf,d_vl_buf,d_ug_buf,d_ud_buf,d_rptr_buf,d_ci_buf,d_vals);
    cudaDeviceSynchronize();
    *d_vals_out=d_vals;
}
