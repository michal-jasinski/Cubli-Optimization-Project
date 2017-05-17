#ifndef __c2_CubliModel_h__
#define __c2_CubliModel_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc2_CubliModelInstanceStruct
#define typedef_SFc2_CubliModelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_isStable;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_CubliModel;
  real_T c2_Theta_0_ht[9];
  real_T c2_m[3];
  real_T (*c2_T)[3];
  real_T (*c2_dot_g)[3];
  real_T (*c2_g)[3];
  real_T (*c2_p_wh)[3];
  real_T (*c2_p_ww)[3];
  real_T (*c2_dot_pwh)[3];
  real_T (*c2_dot_pww)[3];
} SFc2_CubliModelInstanceStruct;

#endif                                 /*typedef_SFc2_CubliModelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_CubliModel_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c2_CubliModel_get_check_sum(mxArray *plhs[]);
extern void c2_CubliModel_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
