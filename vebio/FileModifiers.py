from shutil import copyfile
import os
import re

def write_file_with_replacements(filename, replacements, full_overwrite=False):

    """ Modify an existing file with replacement values.

    This function writes a modified version of an existing file with 
    replacements drawn from a user-provided dictionary.  Candidates
    for replacement are identified by the keyname provided in the dictionary,
    where the new value is stored in ``replacements[keyname]``.  It is written
    to identify variable assignments in the original file only, for example::

        double t_final = 50.0;

    would be the intended target for ``{"t_final": 25.0}`` to halve the variable
    value.  It will preserve the assignment operator for "=", ":", and " " 
    characters, where the last pattern represents any number of blank spaces.
    
    Args:
        filename (str):
            The existing file to modify.

        replacements (str):
            A dictionary of replacement values where the keyname exactly matches
            the value of a variable created with an assignment operator.

        full_overwrite (bool, optional):
            An optional argument to specify that the value associated with 
            ``replacements[keyname]`` should be used to *completely overwrite
            the matching line in the target file*, i.e., it will not attempt 
            to preserve any assignment operators, it will simply overwrite the
            target line with the value from the dictionary with no error checking.
            Defaults to ``False``.

    Returns:
        None

    """


    # Extract the path (if applicable) and the filename
    head, tail = os.path.split(filename)

    # Provide the name of an unchanging backup file to read options from
    backup_filename = os.path.join(head, 'bkup_%s' % (tail))

    try:
        # If it already exists (i.e., we've run this code before) simply 
        # read its contents
        f_read = open(backup_filename, 'r')

    except:
        # If it doesn't exist (i.e., the code is being run for the first time)
        # copy the contents of the original file into a new backup
        # and read the contents from the new backup 
        copyfile(filename, backup_filename)
        f_read = open(backup_filename, 'r')

    # Create a separate pointer for reading and writing
    f_write = open(filename, 'w')

    for line in f_read:

        new_value = None

        # Find the definition of the variable we want to change
        for key, value in replacements.items():
            left = line.split('=')[0]
            split_line = re.split(' |: ', left)
            if key in split_line:
                new_value = value
                del replacements[key]
                break

        if new_value is not None:
            if full_overwrite:
                # If full_overwrite == True, the user must specify a formatted 
                # string that can entirely replace the old line, no error checking
                # on values is done here
                new_line = value

            else:
                # Split on whichever separator/assignment operator is found, and 
                # preserve the left-hand side variable name while adding the
                # new value on the right-hand side
                if '=' in line:
                    lhs = line.split('=')[0]
                    new_line = '%s= %s' % (lhs, str(new_value))

                elif ':' in line:
                    lhs = line.split(':')[0]
                    new_line = '%s: %s' % (lhs, str(new_value))

                else:
                    lhs = line.split()[0]
                    new_line = '%s %s' % (lhs, str(new_value))

                # If the original line ended with a semicolon, add one
                # to the new line as well
                if ';' in line:
                    new_line = '%s;\n' % (new_line)
                else:
                    new_line = '%s\n' % (new_line)

            # Replace the current line with the modified line
            line = new_line

        # Write the (possibly modified) line into the new file    
        f_write.write(line)

    # Close both files when finished
    f_read.close()
    f_write.close()
