#ifndef CMDLINE_H
#define CMDLINE_H
#include "base_types.hh"
#include <string>

// The restricted structure
extern std::string input_structure;
// The parameter file location
extern std::string parameter_file;

/** @brief Where the command line options are stored */
struct args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  const char *verbose_help; /**< @brief Turn on verbose output help description.  */
  const char *input_structure_help; /**< @brief Give restricted structure as input help description.  */
  const char *paramFile_help; /**< @brief Use a separate parameter list */
  const char *noGC_help; /**< @brief Turn off garbage collection and related overhead help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int input_structure_given ;	/**< @brief Whether restricted structure was given.  */
  unsigned int paramFile_given ; /** <@brief whether a parameter file was given */
  unsigned int noGC_given ;	/**< @brief Whether noGC was given.  */

  char **inputs ; /**< @brief unnamed options (options without names) */
  unsigned inputs_num ; /**< @brief unnamed options number */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];


/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,struct args_info *args_info);



/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct args_info *args_info);

/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct args_info *args_info);

#endif