#ifndef __c1_CubliModel_h__
#define __c1_CubliModel_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc1_CubliModelInstanceStruct
#define typedef_SFc1_CubliModelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_CubliModel;
  real_T c1_Theta_0_ht[9];
  real_T c1_m[3];
  real_T c1_alpha;
  real_T c1_beta;
  real_T c1_gamma;
  real_T c1_delta;
  real_T (*c1_g)[3];
  real_T (*c1_T)[3];
  real_T (*c1_p_wh)[3];
  real_T (*c1_p_ww)[3];
} SFc1_CubliModelInstanceStruct;

#endif                                 /*typedef_SFc1_CubliModelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_CubliModel_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_CubliModel_get_check_sum(mxArray *plhs[]);
extern void c1_CubliModel_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
