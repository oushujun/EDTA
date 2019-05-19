/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef OPTION_API_H
#define OPTION_API_H

#include <stdbool.h>
#include "core/range_api.h"
#include "core/str_api.h"
#include "core/str_array_api.h"

/* <GtOptionParser> objects can be used to parse command line options. */
typedef struct GtOptionParser GtOptionParser;
/* <GtOption> objects represent command line options (which are used in
   a <GtOptionParser>).
   Option descriptions are automatically formatted to
   <GT_OPTION_PARSER_TERMINAL_WIDTH>, but it is possible to embed newlines into
   the descriptions to manually affect the formatting. */
typedef struct GtOption GtOption;

enum GtOPrval {
  GT_OPTION_PARSER_OK,           /* Everything went fine. */
  GT_OPTION_PARSER_ERROR,        /* An error occurred during option parsing. */
  GT_OPTION_PARSER_REQUESTS_EXIT /* The option parser requests an exit, because
                                    option -help, -help+, -helpdev, or -version
                                    was used. */
};

/* Possible option parser return values. <GT_OPTION_PARSER_OK> denotes that
   everything went fine, <GT_OPTION_PARSER_ERROR> that an error occurred during
   option parsing, and <GT_OPTION_PARSER_REQUESTS_EXIT> that the option parser
   requests an exit, because option <-help>, <-help+>, <-helpdev> or <-version>
   was used. */
typedef enum GtOPrval GtOPrval;

typedef void (*GtShowVersionFunc)(const char *progname);
typedef int  (*GtShowCommentFunc)(const char *progname, void *data, GtError*);
typedef int  (*GtOptionParserHookFunc)(void *data, GtError*);

/* The default terminal width used in the output of the <GtOptionParser>. */
#define GT_OPTION_PARSER_TERMINAL_WIDTH \
        80

/* Return a new <GtOptionParser> object. The <synopsis> should summarize the
   command line arguments and mandatory arguments in a single line.
   The <one_liner> should describe the program for which the <GtOptionParser> is
   used in a single line and must have an upper case letter at the start and a
   '.' at the end. */
GtOptionParser* gt_option_parser_new(const char *synopsis,
                                     const char *one_liner);
/* Add <option> to <option_parser>. Takes ownership of <option>. */
void            gt_option_parser_add_option(GtOptionParser *option_parser,
                                            GtOption *option);
/* Return the <GtOption> object if an option named <option_string> is present in
   <option_parser>, and <NULL> if no such option exists. */
GtOption*       gt_option_parser_get_option(GtOptionParser *option_parser,
                                            const char *option_string);
/* Refer to manual at the end of <-help> output of <opion_parser>. */
void            gt_option_parser_refer_to_manual(GtOptionParser *option_parser);
/* Set <comment_func> in <option_parser> (<data> is passed along). */
void            gt_option_parser_set_comment_func(GtOptionParser *option_parser,
                                                  GtShowCommentFunc
                                                  comment_func,
                                                  void *data);
/* Set the version function used by <option_parser> to <version_func>.
   This version function takes precedence to the one supplied to
   <gt_option_parser_parse()>. */
void            gt_option_parser_set_version_func(GtOptionParser *option_parser,
                                                  GtShowVersionFunc
                                                  version_func);
/* Set the <mail_address> used in the final "Report bugs to" line of the <-help>
   output. It should be of the form <<bill@microsoft.com>> (email address
   enclosed in one pair of angle brackets). */
void            gt_option_parser_set_mail_address(GtOptionParser*,
                                                  const char *mail_address);
/* Register a <hook_function> with <option_parser>. All registered hook
   functions are called at the end of <gt_option_parser_parse()>.
   This allows one to have a module which registers a bunch of options in the
   option parser and automatically performs necessary postprocessing after the
   option parsing has been done via the hook function. */
void            gt_option_parser_register_hook(GtOptionParser *option_parser,
                                               GtOptionParserHookFunc
                                               hook_function,
                                               void *data);
/* The <minimum> number of additional command line arguments <option_parser>
   must parse in order to succeed. */
void            gt_option_parser_set_min_args(GtOptionParser *option_parser,
                                              unsigned int minimum);
/* The <maximum> number of additional command line arguments <option_parser>
   must parse in order to succeed. */
void            gt_option_parser_set_max_args(GtOptionParser *option_parser,
                                              unsigned int maximum);
/* The <minimum> and <maximum> number of additional command line arguments
   <option_parser> must parse in order to succeed. */
void            gt_option_parser_set_min_max_args(GtOptionParser *option_parser,
                                                  unsigned int minimum,
                                                  unsigned int maximum);
/* Use <option_parser> to parse options given in argument vector <argv> (with
   <argc> many arguments). The number of parsed arguments is stored in
   <parsed_args>. <version_func> is used for the output of option <-version>.
   In case of error, <GT_OPTION_PARSER_ERROR> is returned and <err> is set
   accordingly. */
GtOPrval        gt_option_parser_parse(GtOptionParser *option_parser,
                                       int *parsed_args,
                                       int argc, const char **argv,
                                       GtShowVersionFunc version_func,
                                       GtError *err);
/* Reset all options set in <op> to the default values specified at
   option parser creation time. */
void            gt_option_parser_reset(GtOptionParser *op);
/* Delete <option_parser>. */
void            gt_option_parser_delete(GtOptionParser *option_parser);

/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_bool(const char *option_string,
                                   const char *description,
                                   bool *value, bool default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_double(const char *option_string,
                                     const char *description, double *value,
                                     double default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value>. */
GtOption*       gt_option_new_double_min(const char *option_string,
                                         const char *description, double *value,
                                         double default_value,
                                         double minimum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value> and at most
  the <maximum_value>. */
GtOption*       gt_option_new_double_min_max(const char *option_string,
                                             const char *description,
                                             double *value,
                                             double default_value,
                                             double minimum_value,
                                             double maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at larger or equal than 0.0 and smaller or
  equal than 1.0. */
GtOption*       gt_option_new_probability(const char *option_string,
                                          const char *description,
                                          double *value,
                                          double default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_int(const char *option_string,
                                  const char *description,
                                  int *value, int default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value>. */
GtOption*       gt_option_new_int_min(const char *option_string,
                                      const char *description, int *value,
                                      int default_value, int minimum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at most have the <maximum_value>. */
GtOption*       gt_option_new_int_max(const char *option_string,
                                      const char *description, int *value,
                                      int default_value, int maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value> and at most
  the <maximum_value>. */
GtOption*       gt_option_new_int_min_max(const char *option_string,
                                          const char *description,
                                          int *value, int default_value,
                                          int minimum_value, int maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_uint(const char *option_string,
                                   const char *description,
                                   unsigned int *value,
                                   unsigned int default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value>. */
GtOption*       gt_option_new_uint_min(const char *option_string,
                                       const char *description,
                                       unsigned int *value,
                                       unsigned int default_value,
                                       unsigned int minimum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at most have the <maximum_value>. */
GtOption*       gt_option_new_uint_max(const char *option_string,
                                       const char *description,
                                       unsigned int *value,
                                       unsigned int default_value,
                                       unsigned int maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value> and at most
  the <maximum_value>. */
GtOption*       gt_option_new_uint_min_max(const char *option_string,
                                           const char *description,
                                           unsigned int *value,
                                           unsigned int default_value,
                                           unsigned int minimum_value,
                                           unsigned int maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_word(const char *option_string,
                                   const char *description,
                                   GtWord *value, GtWord default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_uword(const char *option_string,
                                    const char *description,
                                    GtUword *value,
                                    GtUword default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value>. */
GtOption*       gt_option_new_uword_min(const char *option_string,
                                        const char *description,
                                        GtUword *value,
                                        GtUword default_value,
                                        GtUword minimum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The argument to this option must at least have the <minimum_value> and at most
  the <maximum_value>. */
GtOption*       gt_option_new_uword_min_max(const char *option_string,
                                            const char *description,
                                            GtUword *value,
                                            GtUword default_value,
                                            GtUword minimum_value,
                                            GtUword maximum_value);
/* Deprecated. Usage identical to <gt_option_new_word>. */
GtOption*       gt_option_new_long(const char *option_string,
                                   const char *description,
                                   GtWord *value, GtWord default_value);
/* Deprecated. Usage identical to <gt_option_new_uword>. */
GtOption*       gt_option_new_ulong(const char *option_string,
                                    const char *description,
                                    GtUword *value,
                                    GtUword default_value);
/* Deprecated. Usage identical to <gt_option_new_uword_min>. */
GtOption*       gt_option_new_ulong_min(const char *option_string,
                                        const char *description,
                                        GtUword *value,
                                        GtUword default_value,
                                        GtUword minimum_value);
/* Deprecated. Usage identical to <gt_option_new_uword_min_max>. */
GtOption*       gt_option_new_ulong_min_max(const char *option_string,
                                            const char *description,
                                            GtUword *value,
                                            GtUword default_value,
                                            GtUword minimum_value,
                                            GtUword maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
   If <default_value> equals <NULL>, <GT_UNDEF_WORD> will be used as the default
   start and end point of <value>. */
GtOption*       gt_option_new_range(const char *option_string,
                                    const char *description,
                                    GtRange *value, GtRange *default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>.
  The first argument to this option (which will be used as the start) must at
  least have the <minimum_value> and the second argument (which will be used as
  the end) at most the <maximum_value>. */
GtOption*       gt_option_new_range_min_max(const char *option_string,
                                            const char *description,
                                            GtRange *value,
                                            GtRange *default_value,
                                            GtUword minimum_value,
                                            GtUword maximum_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing is stored in <value>. */
GtOption*       gt_option_new_string(const char *option_string,
                                     const char *description,
                                     GtStr *value, const char *default_value);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The result of the option parsing are stored in <value>. */
GtOption*       gt_option_new_string_array(const char *option_string,
                                           const char *description,
                                           GtStrArray *value);
/* Return a <GtOption> with the given <option_string>, <description>, and
   <default_value> which allows only arguments given in the <NULL>-terminated
   <domain> (<default_value> must be an entry of <domain> or <NULL>). */
GtOption*       gt_option_new_choice(const char *option_string,
                                     const char *description, GtStr *value,
                                     const char *default_value,
                                     const char **domain);
/* Return a new <GtOption> with the given <option_string>, <description>. The
   result of the option parsing is stored in the <GtStr> object <filename>.
   <filename> may not be NULL! */
GtOption*       gt_option_new_filename(const char *option_string,
                                       const char *description,
                                       GtStr *filename);
/* Return a new <GtOption> with the given <option_string>, <description>, and
   <default_value>. The results of the option parsing are stored in <value>. */
GtOption*       gt_option_new_filename_array(const char *option_string,
                                             const char *description,
                                             GtStrArray *filename_array);
/* Return a new debug <GtOption> object: <-debug>, "enable debugging output",
   default is <false>. The result of the option parsing is stored in <value> */
GtOption*       gt_option_new_debug(bool *value);
/* Return a new verbose <GtOption> object: <-v>, "be verbose",
   default is <false>. The result of the option parsing is stored in <value> */
GtOption*       gt_option_new_verbose(bool *value);
/* Return a new width <GtOption> object: <-width>, "set output width for FASTA
   sequence printing (0 disables formatting)", default is 0.
   The result of the option parsing is stored in <value> */
GtOption*       gt_option_new_width(GtUword *value);
/* Increase the reference count for <option> and return it. */
GtOption*       gt_option_ref(GtOption *option);
/* Return the name of <option> */
const char*     gt_option_get_name(const GtOption * option);
/* Make <option> mandatory. */
void            gt_option_is_mandatory(GtOption *option);
/* Make it mandatory, that either <option_a> or <option_b> is used. */
void            gt_option_is_mandatory_either(GtOption *option_a,
                                              const GtOption *option_b);
/* Make it mandatory, that one of the options <option_a>, <option_b>, or
   <option_c> is used. */
void            gt_option_is_mandatory_either_3(GtOption *option_a,
                                                const GtOption *option_b,
                                                const GtOption *option_c);
/* Make it mandatory, that one of the options <option_a>, <option_b>, <option_c>
   or <option_d> is used. */
void            gt_option_is_mandatory_either_4(GtOption *option_a,
                                                const GtOption *option_b,
                                                const GtOption *option_c,
                                                const GtOption *option_d);
/* Set that <option> is only shown in the output of <-help+>. */
void            gt_option_is_extended_option(GtOption *option);
/* Set that <option> is only shown in the output of <-helpdev>. */
void            gt_option_is_development_option(GtOption *option);
/* Make <option_a> imply <option_b>. */
void            gt_option_imply(GtOption *option_a, const GtOption *option_b);
/* Make <option_a> imply either <option_b> or <option_c> */
void            gt_option_imply_either_2(GtOption *option_a,
                                         const GtOption *option_b,
                                         const GtOption *option_c);
/* Make <option_a> imply either <option_b>, <option_c> or <option_d> */
void            gt_option_imply_either_3(GtOption *option_a,
                                         const GtOption *option_b,
                                         const GtOption *option_c,
                                         const GtOption *option_d);
/* Set that the options <option_a> and <option_b> exclude each other. */
void            gt_option_exclude(GtOption *option_a, GtOption *option_b);
/* Hide the default value of <option> in <-help> output. */
void            gt_option_hide_default(GtOption *option);
/* Set that the argument to <option> is optional */
void            gt_option_argument_is_optional(GtOption *option);
/* Return <true> if <option> was set, <false> otherwise. */
bool            gt_option_is_set(const GtOption *option);
/* Delete <option>. */
void            gt_option_delete(GtOption*);

/* Parse the argument to option -memlimit. Could be made into
   a special parser, but I do not know how. SK. 2011-09-19 */

int gt_option_parse_spacespec(GtUword *maximumspace,
                              const char *optname,
                              const GtStr *memlimit,
                              GtError *err);

#endif
