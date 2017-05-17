#ifndef __c3_CubliModel_h__
#define __c3_CubliModel_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc3_CubliModelInstanceStruct
#define typedef_SFc3_CubliModelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c3_sfEvent;
  boolean_T c3_isStable;
  boolean_T c3_doneDoubleBufferReInit;
  uint8_T c3_is_active_c3_CubliModel;
  real_T c3_Theta_0_ht[9];
  real_T (*c3_p_wh)[3];
  real_T (*c3_PIK_dot)[4];
  real_T (*c3_p_ww)[3];
  real_T (*c3_PIK)[4];
} SFc3_CubliModelInstanceStruct;

#endif                                 /*typedef_SFc3_CubliModelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c3_CubliModel_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c3_CubliModel_get_check_sum(mxArray *plhs[]);
extern void c3_CubliModel_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
