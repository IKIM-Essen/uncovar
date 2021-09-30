
<details><summary>Expected commit message structure</summary>

Commit messages should be structured as follows:

```
<type>[optional scope]: <description>

[optional body]

[optional footer]
```

## Type
Must be one of the following:

- **build**: Changes that affect the build system or external dependencies (example scopes: gulp, broccoli, npm)
- **ci**: Changes to our CI configuration files and scripts (example scopes: Travis, Circle, BrowserStack, SauceLabs)
- **docs**: Documentation only changes
- **feat**: A new feature
- **fix**: A bug fix
- **perf**: A code change that improves performance
- **refactor**: A code change that neither fixes a bug nor adds a feature
- **style**: Changes that do not affect the meaning of the code (white-space, formatting, missing semi-colons, etc)
- **test**: Adding missing tests or correcting existing tests

## Examples

Commit message with description and breaking change in body
```
feat: allow provided config object to extend other configs

BREAKING CHANGE: `extends` key in config file is now used for extending other config files
```

Commit message with no body
```
docs: correct spelling of CHANGELOG
````

Commit message with scope
```
feat(lang): added polish language
```

Commit message for a fix using an (optional) issue number.
```
fix: minor typos in code

see the issue for details on the typos fixed

fixes issue #12
```

</p>
</details>
