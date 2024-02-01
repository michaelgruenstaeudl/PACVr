import os
import sys
import re


def main():
    r_directory = "../../R/"
    current_timestamp = sys.argv[1]
    version_comment = f"#version=\"{current_timestamp}\"\n"

    for filename in os.listdir(r_directory):
        file_path = os.path.join(r_directory, filename)

        if os.path.isfile(file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()

            if len(lines) >= 4 and re.match("^#version", lines[3]):
                lines[3] = version_comment
                with open(file_path, 'w') as file:
                    file.writelines(lines)
                    print(f"Version comment updated to current UTC: {filename}")


if __name__ == "__main__":
    main()
