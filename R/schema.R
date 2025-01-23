#' This function will print out infromation about a given schema Key
#'
#' @param key The schema Key to provide information about to the user.
#'
#' @return
#'
#' @export
._get_schema_def = function(key=NULL){
    return(schema %>% filter(Key == key))
}

# .validate_schema = function(){
#     """
#         Given a schema definition and the column of a dataframe, verify that the column satisfies the schema.
#         If there are any errors, it adds them to self.errors

#         :param pandas.Series column: Column of a dataframe to validate
#         :param str column_name: Name of the column in the dataframe
#         :param str df_name: Name of the dataframe
#         :param dict column_def: schema definition for this specific column,
#         e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]
#         :param str default_error_message_suffix: default error message suffix to be added to errors found here

#         :rtype None
#     """
# }

# ._validate_obs_column(column_name, column_def){
#     ## Validate the column
#     if(column_def$Type == "str")
#         assert(all(is.character(column)), "Column must be boolean.")

#     if(column_def$Type == "bool")
#         assert(all(column %in% c(TRUE, FALSE)), "Column must be boolean.")
# }

#     def _validate_column(
#         self, column: pd.Series, column_name: str, df_name: str, column_def: dict, default_error_message_suffix=None
#     ):
#         """
#         Given a schema definition and the column of a dataframe, verify that the column satisfies the schema.
#         If there are any errors, it adds them to self.errors

#         :param pandas.Series column: Column of a dataframe to validate
#         :param str column_name: Name of the column in the dataframe
#         :param str df_name: Name of the dataframe
#         :param dict column_def: schema definition for this specific column,
#         e.g. schema_def["obs"]["columns"]["cell_type_ontology_term_id"]
#         :param str default_error_message_suffix: default error message suffix to be added to errors found here

#         :rtype None
#         """

#         # error_original_count will count the number of error messages prior to validating the column, this
#         # will be useful in case there's an error prefix to be added to errors found here
#         error_original_count = len(self.errors)

#         if column_def.get("unique") and column.nunique() != len(column):
#             self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' is not unique.")

#         if column_def.get("type") == "bool" and column.dtype != bool:
#             self.errors.append(
#                 f"Column '{column_name}' in dataframe '{df_name}' must be boolean, not '{column.dtype.name}'."
#             )

#         if column_def.get("type") == "categorical":
#             if column.dtype.name != "category":
#                 self.errors.append(
#                     f"Column '{column_name}' in dataframe '{df_name}' must be categorical, not {column.dtype.name}."
#                 )
#             else:
#                 if column_def.get("subtype") == "str":
#                     if column.dtype.categories.dtype != "object" and column.dtype.categories.dtype != "string":
#                         self.errors.append(
#                             f"Column '{column_name}' in dataframe '{df_name}' must be object or string, not"
#                             f" {column.dtype.categories.dtype}."
#                         )
#                     else:
#                         if any(len(cat.strip()) == 0 for cat in column.dtype.categories):
#                             self.errors.append(
#                                 f"Column '{column_name}' in dataframe '{df_name}' must not contain empty values."
#                             )

#                 # check for null values--skip on column defs with enums, since it will already be part of that check
#                 if not column_def.get("enum") and column.isnull().any():
#                     self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' must not contain NaN values.")

#         if column_def.get("type") == "feature_is_filtered":
#             self._validate_column_feature_is_filtered(column, column_name, df_name)

#         if column_def.get("type") == "genetic_ancestry_value":
#             self._validate_individual_genetic_ancestry_value(column, column_name)

#         if "enum" in column_def:
#             bad_enums = [v for v in column.drop_duplicates() if v not in column_def["enum"]]
#             if bad_enums:
#                 self.errors.append(
#                     f"Column '{column_name}' in dataframe '{df_name}' contains invalid values "
#                     f"'{bad_enums}'. Values must be one of {column_def['enum']}"
#                 )

#         if column_def.get("type") == "feature_id":
#             # Validates each id
#             for feature_id in column:
#                 self._validate_feature_id(feature_id, df_name)

#         if column_def.get("type") == "curie":
#             # Check for NaN values
#             if column.isnull().any():
#                 self.errors.append(f"Column '{column_name}' in dataframe '{df_name}' must not contain NaN values.")
#                 return

#             if "curie_constraints" in column_def:
#                 for term_str in column.drop_duplicates():
#                     self._validate_curie_str(term_str, column_name, column_def["curie_constraints"])

#         # Add error suffix to errors found here
#         error_message_suffix = column_def.get("error_message_suffix", default_error_message_suffix)
#         if error_message_suffix:
#             error_total_count = len(self.errors)
#             for i in range(error_original_count, error_total_count):
#                 self.errors[i] = self.errors[i] + " " + error_message_suffix